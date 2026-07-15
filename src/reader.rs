use crate::cli::PhredOffset;
use crate::metrics::{
    MAX_TILES, MAX_TRACKED_UNIQUE_SEQUENCES, NamedSequence, OVERREP_FAIL_PERCENT,
    OVERREP_WARN_PERCENT, PositionGroup, high_is_bad, low_is_bad, make_base_groups,
    possible_source,
};
use crate::report::{
    AdapterContentData, AdapterSeries, BaseContentPoint, BaseQualityPoint, BasicStatistics,
    CountPoint, DuplicationData, DuplicationPoint, GcContentData, GcPoint, Module,
    OverrepresentedSequence, PositionValue, Report, SequenceQualityData, Status, TileQualityData,
};
use crossbeam_channel::{Receiver, Sender, TrySendError, bounded};
use indicatif::{ProgressBar, ProgressStyle};
use needletail::{FastxReader, parse_fastx_reader, parse_fastx_stdin};
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::Read;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::thread;
use std::time::{SystemTime, UNIX_EPOCH};

const QUALITY_MIN_ASCII: usize = 33;
const QUALITY_MAX_ASCII: usize = 126;
const QUALITY_BINS: usize = QUALITY_MAX_ASCII - QUALITY_MIN_ASCII + 1;
const RECORDS_PER_BATCH: usize = 4_096;

#[derive(Clone)]
pub struct AnalysisConfig {
    pub no_group: bool,
    pub nofilter: bool,
    pub dup_length: usize,
    pub min_length: usize,
    pub phred_offset: PhredOffset,
    pub adapters: Vec<NamedSequence>,
    pub contaminants: Vec<NamedSequence>,
}

#[derive(Clone)]
struct PositionAccumulator {
    quality_counts: [u64; QUALITY_BINS],
    quality_sum_ascii: u64,
    quality_count: u64,
    bases: [u64; 5],
}

impl Default for PositionAccumulator {
    fn default() -> Self {
        Self {
            quality_counts: [0; QUALITY_BINS],
            quality_sum_ascii: 0,
            quality_count: 0,
            bases: [0; 5],
        }
    }
}

#[derive(Clone, Copy, Default)]
struct SumCount {
    sum: u64,
    count: u64,
}

struct CountingReader<R> {
    inner: R,
    bytes_read: Arc<AtomicU64>,
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let count = self.inner.read(buffer)?;
        self.bytes_read.fetch_add(count as u64, Ordering::Relaxed);
        Ok(count)
    }
}

struct OwnedRecord {
    id: Vec<u8>,
    sequence: Vec<u8>,
    quality: Option<Vec<u8>>,
}

struct WorkBatch {
    index: usize,
    records: Vec<OwnedRecord>,
}

struct OrderedRecord {
    duplicate_prefix: Vec<u8>,
    tile: Option<(Option<u32>, Vec<u8>)>,
}

struct BatchResult {
    index: usize,
    accumulator: Accumulator,
    ordered_records: Vec<OrderedRecord>,
}

struct BatchFailure {
    index: usize,
    error: String,
}

type BatchOutcome = Result<BatchResult, BatchFailure>;

struct Accumulator {
    filename: String,
    total_sequences: u64,
    analyzed_sequences: u64,
    total_bases: u64,
    total_gc: u64,
    poor_quality_sequences: u64,
    min_length: usize,
    max_length: usize,
    length_counts: BTreeMap<usize, u64>,
    positions: Vec<PositionAccumulator>,
    min_quality_ascii: u8,
    has_quality: bool,
    sequence_quality_ascii: [u64; 128],
    gc_distribution: [f64; 101],
    gc_claim_cache: HashMap<usize, [u32; 101]>,
    sequence_counts: HashMap<Vec<u8>, u64>,
    duplication_frozen: bool,
    count_at_unique_limit: u64,
    tile_counts: HashMap<u32, Vec<SumCount>>,
    tile_observations: u64,
    tile_supported: bool,
    adapter_starts: Vec<Vec<u64>>,
}

impl Accumulator {
    fn new(filename: String, adapter_count: usize) -> Self {
        Self::with_sequence_capacity(filename, adapter_count, MAX_TRACKED_UNIQUE_SEQUENCES)
    }

    fn batch(filename: String, adapter_count: usize) -> Self {
        Self::with_sequence_capacity(filename, adapter_count, 0)
    }

    fn with_sequence_capacity(
        filename: String,
        adapter_count: usize,
        sequence_capacity: usize,
    ) -> Self {
        Self {
            filename,
            total_sequences: 0,
            analyzed_sequences: 0,
            total_bases: 0,
            total_gc: 0,
            poor_quality_sequences: 0,
            min_length: usize::MAX,
            max_length: 0,
            length_counts: BTreeMap::new(),
            positions: Vec::new(),
            min_quality_ascii: u8::MAX,
            has_quality: false,
            sequence_quality_ascii: [0; 128],
            gc_distribution: [0.0; 101],
            gc_claim_cache: HashMap::new(),
            sequence_counts: HashMap::with_capacity(sequence_capacity),
            duplication_frozen: false,
            count_at_unique_limit: 0,
            tile_counts: HashMap::new(),
            tile_observations: 0,
            tile_supported: true,
            adapter_starts: vec![Vec::new(); adapter_count],
        }
    }

    fn process_record(
        &mut self,
        id: &[u8],
        sequence: &[u8],
        quality: Option<&[u8]>,
        config: &AnalysisConfig,
    ) -> Result<(), String> {
        if !self.process_reducible_record(id, sequence, quality, config)? {
            return Ok(());
        }
        self.add_duplication_observation(sequence, config.dup_length);
        if let Some(quality) = quality.filter(|quality| !quality.is_empty()) {
            self.add_tile_observation(id, quality);
        }
        Ok(())
    }

    /// Updates statistics that can be merged exactly from independently analyzed batches.
    fn process_reducible_record(
        &mut self,
        id: &[u8],
        sequence: &[u8],
        quality: Option<&[u8]>,
        config: &AnalysisConfig,
    ) -> Result<bool, String> {
        self.total_sequences += 1;
        self.total_bases += sequence.len() as u64;
        self.total_gc += sequence
            .iter()
            .filter(|base| matches!(base.to_ascii_uppercase(), b'G' | b'C'))
            .count() as u64;
        self.min_length = self.min_length.min(sequence.len());
        self.max_length = self.max_length.max(sequence.len());
        *self.length_counts.entry(sequence.len()).or_default() += 1;

        let poor_quality = is_casava_filtered(id);
        if poor_quality {
            self.poor_quality_sequences += 1;
        }
        if poor_quality && !config.nofilter {
            return Ok(false);
        }
        self.analyzed_sequences += 1;

        if self.positions.len() < sequence.len() {
            self.positions
                .resize_with(sequence.len(), PositionAccumulator::default);
        }

        let mut gc_count = 0usize;
        for (index, raw_base) in sequence.iter().copied().enumerate() {
            let base = raw_base.to_ascii_uppercase();
            let slot = match base {
                b'A' => 0,
                b'C' => {
                    gc_count += 1;
                    1
                }
                b'G' => {
                    gc_count += 1;
                    2
                }
                b'T' | b'U' => 3,
                _ => 4,
            };
            self.positions[index].bases[slot] += 1;
        }
        self.add_gc_observation(sequence, gc_count);
        self.add_adapter_observations(sequence, &config.adapters);

        if let Some(quality) = quality {
            if quality.len() != sequence.len() {
                return Err(format!(
                    "{}: sequence and quality lengths differ for record {}",
                    self.filename,
                    String::from_utf8_lossy(id)
                ));
            }
            if !quality.is_empty() {
                self.has_quality = true;
                let mut quality_sum = 0u64;
                for (index, value) in quality.iter().copied().enumerate() {
                    if !(QUALITY_MIN_ASCII as u8..=QUALITY_MAX_ASCII as u8).contains(&value) {
                        return Err(format!(
                            "{}: unsupported quality byte {value} in record {}",
                            self.filename,
                            String::from_utf8_lossy(id)
                        ));
                    }
                    self.min_quality_ascii = self.min_quality_ascii.min(value);
                    quality_sum += value as u64;
                    let position = &mut self.positions[index];
                    position.quality_counts[value as usize - QUALITY_MIN_ASCII] += 1;
                    position.quality_sum_ascii += value as u64;
                    position.quality_count += 1;
                }
                let average_ascii = (quality_sum / quality.len() as u64) as usize;
                self.sequence_quality_ascii[average_ascii] += 1;
            }
        }
        Ok(true)
    }

    fn add_gc_observation(&mut self, sequence: &[u8], full_gc_count: usize) {
        if sequence.is_empty() {
            return;
        }
        let model_length = if sequence.len() > 1000 {
            (sequence.len() / 1000) * 1000
        } else if sequence.len() > 100 {
            (sequence.len() / 100) * 100
        } else {
            sequence.len()
        };
        let gc_count = if model_length == sequence.len() {
            full_gc_count
        } else {
            sequence[..model_length]
                .iter()
                .filter(|base| matches!(base.to_ascii_uppercase(), b'G' | b'C'))
                .count()
        };
        let claims = self
            .gc_claim_cache
            .entry(model_length)
            .or_insert_with(|| gc_claim_counts(model_length));
        let (low, high) = gc_percentage_range(gc_count, model_length);
        for (percentage, claim_count) in claims.iter().enumerate().take(high + 1).skip(low) {
            self.gc_distribution[percentage] += 1.0 / *claim_count as f64;
        }
    }

    fn add_duplication_observation(&mut self, sequence: &[u8], duplication_length: usize) {
        let length = sequence.len().min(duplication_length);
        let prefix = &sequence[..length];
        if let Some(count) = self.sequence_counts.get_mut(prefix) {
            *count += 1;
            if !self.duplication_frozen {
                self.count_at_unique_limit = self.analyzed_sequences;
            }
            return;
        }
        if prefix.iter().any(u8::is_ascii_lowercase) {
            self.add_duplication_prefix(prefix.to_ascii_uppercase());
        } else {
            self.add_duplication_prefix(prefix.to_vec());
        }
    }

    fn add_duplication_prefix(&mut self, prefix: Vec<u8>) {
        self.add_duplication_prefix_at(prefix, self.analyzed_sequences);
    }

    fn add_duplication_prefix_at(&mut self, prefix: Vec<u8>, analyzed_position: u64) {
        if let Some(count) = self.sequence_counts.get_mut(prefix.as_slice()) {
            *count += 1;
            if !self.duplication_frozen {
                self.count_at_unique_limit = analyzed_position;
            }
            return;
        }
        if !self.duplication_frozen {
            self.sequence_counts.insert(prefix, 1);
            self.count_at_unique_limit = analyzed_position;
            if self.sequence_counts.len() == MAX_TRACKED_UNIQUE_SEQUENCES {
                self.duplication_frozen = true;
            }
        }
    }

    fn add_tile_observation(&mut self, id: &[u8], quality: &[u8]) {
        self.add_tile_observation_for_tile(parse_tile(id), quality);
    }

    fn add_tile_observation_for_tile(&mut self, tile: Option<u32>, quality: &[u8]) {
        if !self.tile_supported {
            return;
        }
        self.tile_observations += 1;
        if self.tile_observations > 10_000 && !self.tile_observations.is_multiple_of(10) {
            return;
        }
        let Some(tile) = tile else {
            self.tile_supported = false;
            self.tile_counts.clear();
            return;
        };
        if !self.tile_counts.contains_key(&tile) && self.tile_counts.len() >= MAX_TILES {
            self.tile_supported = false;
            self.tile_counts.clear();
            return;
        }
        let counts = self.tile_counts.entry(tile).or_default();
        if counts.len() < quality.len() {
            counts.resize(quality.len(), SumCount::default());
        }
        for (position, value) in counts.iter_mut().zip(quality) {
            position.sum += *value as u64;
            position.count += 1;
        }
    }

    fn add_adapter_observations(&mut self, sequence: &[u8], adapters: &[NamedSequence]) {
        for (adapter_index, adapter) in adapters.iter().enumerate() {
            if let Some(position) = find_subslice_case_insensitive(sequence, &adapter.sequence) {
                let starts = &mut self.adapter_starts[adapter_index];
                if starts.len() <= position {
                    starts.resize(position + 1, 0);
                }
                starts[position] += 1;
            }
        }
    }

    fn merge_reducible(&mut self, partial: Accumulator) {
        self.total_sequences += partial.total_sequences;
        self.analyzed_sequences += partial.analyzed_sequences;
        self.total_bases += partial.total_bases;
        self.total_gc += partial.total_gc;
        self.poor_quality_sequences += partial.poor_quality_sequences;
        self.min_length = self.min_length.min(partial.min_length);
        self.max_length = self.max_length.max(partial.max_length);
        for (length, count) in partial.length_counts {
            *self.length_counts.entry(length).or_default() += count;
        }
        if self.positions.len() < partial.positions.len() {
            self.positions
                .resize_with(partial.positions.len(), PositionAccumulator::default);
        }
        for (target, source) in self.positions.iter_mut().zip(partial.positions) {
            for (target_count, source_count) in
                target.quality_counts.iter_mut().zip(source.quality_counts)
            {
                *target_count += source_count;
            }
            target.quality_sum_ascii += source.quality_sum_ascii;
            target.quality_count += source.quality_count;
            for (target_count, source_count) in target.bases.iter_mut().zip(source.bases) {
                *target_count += source_count;
            }
        }
        self.min_quality_ascii = self.min_quality_ascii.min(partial.min_quality_ascii);
        self.has_quality |= partial.has_quality;
        for (target, source) in self
            .sequence_quality_ascii
            .iter_mut()
            .zip(partial.sequence_quality_ascii)
        {
            *target += source;
        }
        for (target, source) in self.gc_distribution.iter_mut().zip(partial.gc_distribution) {
            *target += source;
        }
        for (target, source) in self.adapter_starts.iter_mut().zip(partial.adapter_starts) {
            if target.len() < source.len() {
                target.resize(source.len(), 0);
            }
            for (target_count, source_count) in target.iter_mut().zip(source) {
                *target_count += source_count;
            }
        }
    }
}

#[cfg(test)]
pub fn analyze_input(input: &str, config: &AnalysisConfig) -> Result<Report, String> {
    analyze_input_with_workers(input, config, 1, false)
}

pub fn analyze_input_with_workers(
    input: &str,
    config: &AnalysisConfig,
    workers: usize,
    show_progress: bool,
) -> Result<Report, String> {
    let filename = input_filename(input);
    let mut accumulator = Accumulator::new(filename, config.adapters.len());
    let progress = make_progress(input, show_progress);
    let byte_counter = Arc::new(AtomicU64::new(0));
    if input == "-" || input.starts_with("stdin:") {
        let mut reader = parse_fastx_stdin().map_err(|error| format!("stdin: {error}"))?;
        consume_with_workers(
            &mut *reader,
            &mut accumulator,
            config,
            workers,
            progress.as_ref(),
            None,
        )?;
    } else {
        let file = File::open(input).map_err(|error| format!("failed to open {input}: {error}"))?;
        let source = CountingReader {
            inner: file,
            bytes_read: Arc::clone(&byte_counter),
        };
        let mut reader = parse_fastx_reader(source)
            .map_err(|error| format!("failed to open {input}: {error}"))?;
        consume_with_workers(
            &mut *reader,
            &mut accumulator,
            config,
            workers,
            progress.as_ref(),
            Some(byte_counter.as_ref()),
        )?;
    }
    if accumulator.total_sequences == 0 {
        if let Some(progress) = &progress {
            progress.abandon_with_message("no sequence records found");
        }
        return Err(format!(
            "{} contains no sequence records",
            accumulator.filename
        ));
    }
    let report = finalize(accumulator, config)?;
    if let Some(progress) = &progress {
        progress.finish_with_message("analysis complete");
    }
    Ok(report)
}

fn make_progress(input: &str, enabled: bool) -> Option<ProgressBar> {
    if !enabled {
        return None;
    }
    let total = std::fs::metadata(input).ok().map(|metadata| metadata.len());
    let progress = total
        .filter(|size| *size > 0)
        .map_or_else(ProgressBar::new_spinner, ProgressBar::new);
    let template = if total.is_some() {
        "{spinner:.cyan} {msg:28} [{bar:40.cyan/blue}] {bytes}/{total_bytes}"
    } else {
        "{spinner:.cyan} {msg} {pos} reads"
    };
    let style = ProgressStyle::with_template(template)
        .unwrap_or_else(|_| ProgressStyle::default_spinner())
        .tick_strings(&["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"])
        .progress_chars("=>-");
    progress.set_style(style);
    progress.set_message(format!("analyzing {}", input_filename(input)));
    progress.enable_steady_tick(std::time::Duration::from_millis(80));
    Some(progress)
}

fn update_progress(progress: Option<&ProgressBar>, byte_counter: Option<&AtomicU64>, records: u64) {
    let Some(progress) = progress else {
        return;
    };
    progress.set_position(byte_counter.map_or(records, |counter| counter.load(Ordering::Relaxed)));
    progress.set_message(format!("analyzing ({records} reads)"));
}

fn consume_with_workers(
    reader: &mut dyn FastxReader,
    accumulator: &mut Accumulator,
    config: &AnalysisConfig,
    workers: usize,
    progress: Option<&ProgressBar>,
    byte_counter: Option<&AtomicU64>,
) -> Result<(), String> {
    if workers <= 1 {
        return consume_reader_sequential(reader, accumulator, config, progress, byte_counter);
    }
    consume_reader_parallel(reader, accumulator, config, workers, progress, byte_counter)
}

fn consume_reader_sequential(
    reader: &mut dyn FastxReader,
    accumulator: &mut Accumulator,
    config: &AnalysisConfig,
    progress: Option<&ProgressBar>,
    byte_counter: Option<&AtomicU64>,
) -> Result<(), String> {
    while let Some(record) = reader.next() {
        let record = record.map_err(|error| format!("{}: {error}", accumulator.filename))?;
        let sequence = record.seq();
        accumulator.process_record(record.id(), sequence.as_ref(), record.qual(), config)?;
        if accumulator
            .total_sequences
            .is_multiple_of(RECORDS_PER_BATCH as u64)
        {
            update_progress(progress, byte_counter, accumulator.total_sequences);
        }
    }
    update_progress(progress, byte_counter, accumulator.total_sequences);
    Ok(())
}

fn consume_reader_parallel(
    reader: &mut dyn FastxReader,
    accumulator: &mut Accumulator,
    config: &AnalysisConfig,
    workers: usize,
    progress: Option<&ProgressBar>,
    byte_counter: Option<&AtomicU64>,
) -> Result<(), String> {
    let queue_depth = workers.saturating_mul(2).clamp(2, 32);
    let (work_tx, work_rx) = bounded::<WorkBatch>(queue_depth);
    let (result_tx, result_rx) = bounded::<BatchOutcome>(queue_depth);
    let filename = accumulator.filename.clone();
    let mut pending = BTreeMap::new();
    let mut next_batch = 0usize;
    let mut sent = 0usize;
    let mut received = 0usize;
    let mut failure = None;
    let mut parse_failure = None;

    thread::scope(|scope| {
        for _ in 0..workers {
            let work_rx = work_rx.clone();
            let result_tx = result_tx.clone();
            scope.spawn(|| worker_loop(work_rx, result_tx, config, &filename));
        }
        drop(result_tx);

        let mut batch_index = 0usize;
        'read: loop {
            let mut records = Vec::with_capacity(RECORDS_PER_BATCH);
            let mut reached_eof = false;
            while records.len() < RECORDS_PER_BATCH {
                let Some(record) = reader.next() else {
                    reached_eof = true;
                    break;
                };
                match record {
                    Ok(record) => {
                        let sequence = record.seq();
                        records.push(OwnedRecord {
                            id: record.id().to_vec(),
                            sequence: sequence.into_owned(),
                            quality: record.qual().map(ToOwned::to_owned),
                        });
                    }
                    Err(error) => {
                        parse_failure = Some(format!("{}: {error}", accumulator.filename));
                        break;
                    }
                }
            }
            if !records.is_empty() {
                if let Err(error) = send_batch(
                    WorkBatch {
                        index: batch_index,
                        records,
                    },
                    &work_tx,
                    &result_rx,
                    accumulator,
                    &mut pending,
                    &mut next_batch,
                    &mut received,
                    &mut failure,
                ) {
                    failure.get_or_insert(error);
                    break 'read;
                }
                sent += 1;
                batch_index += 1;
                update_progress(
                    progress,
                    byte_counter,
                    (batch_index * RECORDS_PER_BATCH) as u64,
                );
            }
            if parse_failure.is_some() || reached_eof {
                break;
            }
        }
        drop(work_tx);
        while received < sent {
            if let Err(error) = receive_one_result(
                &result_rx,
                accumulator,
                &mut pending,
                &mut next_batch,
                &mut received,
                &mut failure,
            ) {
                failure.get_or_insert(error);
                break;
            }
        }
    });

    update_progress(progress, byte_counter, accumulator.total_sequences);
    if let Some(error) = parse_failure.or(failure) {
        Err(error)
    } else {
        Ok(())
    }
}

fn worker_loop(
    work_rx: Receiver<WorkBatch>,
    result_tx: Sender<BatchOutcome>,
    config: &AnalysisConfig,
    filename: &str,
) {
    for batch in work_rx {
        let index = batch.index;
        let outcome =
            process_batch(batch, config, filename).map_err(|error| BatchFailure { index, error });
        if result_tx.send(outcome).is_err() {
            break;
        }
    }
}

fn process_batch(
    batch: WorkBatch,
    config: &AnalysisConfig,
    filename: &str,
) -> Result<BatchResult, String> {
    let mut accumulator = Accumulator::batch(filename.to_string(), config.adapters.len());
    let mut ordered_records = Vec::with_capacity(batch.records.len());
    for OwnedRecord {
        id,
        sequence,
        quality,
    } in batch.records
    {
        if !accumulator.process_reducible_record(&id, &sequence, quality.as_deref(), config)? {
            continue;
        }
        let duplicate_prefix =
            sequence[..sequence.len().min(config.dup_length)].to_ascii_uppercase();
        let tile = quality
            .filter(|quality| !quality.is_empty())
            .map(|quality| (parse_tile(&id), quality));
        ordered_records.push(OrderedRecord {
            duplicate_prefix,
            tile,
        });
    }
    Ok(BatchResult {
        index: batch.index,
        accumulator,
        ordered_records,
    })
}

#[allow(clippy::too_many_arguments)]
fn send_batch(
    mut batch: WorkBatch,
    work_tx: &Sender<WorkBatch>,
    result_rx: &Receiver<BatchOutcome>,
    accumulator: &mut Accumulator,
    pending: &mut BTreeMap<usize, BatchResult>,
    next_batch: &mut usize,
    received: &mut usize,
    failure: &mut Option<String>,
) -> Result<(), String> {
    loop {
        match work_tx.try_send(batch) {
            Ok(()) => return Ok(()),
            Err(TrySendError::Full(returned)) => {
                batch = returned;
                receive_one_result(
                    result_rx,
                    accumulator,
                    pending,
                    next_batch,
                    received,
                    failure,
                )?;
            }
            Err(TrySendError::Disconnected(_)) => {
                return Err("analysis workers stopped unexpectedly".into());
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn receive_one_result(
    result_rx: &Receiver<BatchOutcome>,
    accumulator: &mut Accumulator,
    pending: &mut BTreeMap<usize, BatchResult>,
    next_batch: &mut usize,
    received: &mut usize,
    failure: &mut Option<String>,
) -> Result<(), String> {
    let outcome = result_rx
        .recv()
        .map_err(|_| "analysis workers stopped unexpectedly".to_string())?;
    *received += 1;
    accept_batch_result(outcome, accumulator, pending, next_batch, failure);
    while let Ok(outcome) = result_rx.try_recv() {
        *received += 1;
        accept_batch_result(outcome, accumulator, pending, next_batch, failure);
    }
    Ok(())
}

fn accept_batch_result(
    outcome: BatchOutcome,
    accumulator: &mut Accumulator,
    pending: &mut BTreeMap<usize, BatchResult>,
    next_batch: &mut usize,
    failure: &mut Option<String>,
) {
    match outcome {
        Ok(result) => {
            pending.insert(result.index, result);
        }
        Err(error) => {
            failure.get_or_insert(format!("batch {}: {}", error.index, error.error));
        }
    }
    if failure.is_some() {
        return;
    }
    while let Some(result) = pending.remove(next_batch) {
        merge_batch_result(accumulator, result);
        *next_batch += 1;
    }
}

fn merge_batch_result(accumulator: &mut Accumulator, result: BatchResult) {
    let analyzed_before = accumulator.analyzed_sequences;
    accumulator.merge_reducible(result.accumulator);
    for (index, record) in result.ordered_records.into_iter().enumerate() {
        accumulator
            .add_duplication_prefix_at(record.duplicate_prefix, analyzed_before + index as u64 + 1);
        if let Some((tile, quality)) = record.tile {
            accumulator.add_tile_observation_for_tile(tile, &quality);
        }
    }
}

fn finalize(acc: Accumulator, config: &AnalysisConfig) -> Result<Report, String> {
    let (quality_offset, encoding) = quality_encoding(&acc, config.phred_offset)?;
    let groups = make_base_groups(acc.max_length, config.no_group, config.min_length);
    let per_base_quality = finish_per_base_quality(&acc, &groups, quality_offset);
    let per_base_quality_status = quality_status(&per_base_quality, acc.has_quality);
    let tile_data = finish_tile_quality(&acc, &groups);
    let tile_status = if !acc.has_quality || !acc.tile_supported || acc.tile_counts.is_empty() {
        Status::Skip
    } else {
        high_is_bad(tile_data.max_deviation, 5.0, 10.0)
    };
    let sequence_quality = finish_sequence_quality(&acc, quality_offset);
    let sequence_quality_status = if !acc.has_quality {
        Status::Skip
    } else {
        low_is_bad(sequence_quality.most_frequent_score as f64, 28.0, 21.0)
    };
    let base_content = finish_base_content(&acc, &groups);
    let base_content_status = base_content
        .iter()
        .map(|point| (point.g - point.c).abs().max((point.a - point.t).abs()))
        .fold(0.0, f64::max);
    let base_content_status = high_is_bad(base_content_status, 10.0, 20.0);
    let gc_content = finish_gc_content(&acc.gc_distribution);
    let gc_status = high_is_bad(gc_content.deviation_percent, 15.0, 30.0);
    let n_content = finish_n_content(&acc, &groups);
    let n_status = high_is_bad(
        n_content
            .iter()
            .map(|point| point.value)
            .fold(0.0, f64::max),
        5.0,
        20.0,
    );
    let lengths = finish_lengths(&acc.length_counts);
    let length_status = if acc.length_counts.contains_key(&0) {
        Status::Fail
    } else if acc.length_counts.len() > 1 {
        Status::Warn
    } else {
        Status::Pass
    };
    let duplication = finish_duplication(
        &acc.sequence_counts,
        acc.count_at_unique_limit,
        acc.analyzed_sequences,
    );
    let duplication_status = low_is_bad(duplication.total_deduplicated_percent, 70.0, 50.0);
    let overrepresented = finish_overrepresented(
        &acc.sequence_counts,
        acc.analyzed_sequences,
        &config.contaminants,
    );
    let overrepresented_status = overrepresented.first().map_or(Status::Pass, |top| {
        if top.percentage > OVERREP_FAIL_PERCENT {
            Status::Fail
        } else {
            Status::Warn
        }
    });
    let adapters = finish_adapters(&acc, config, &groups);
    let adapter_status = high_is_bad(adapters.max_content, 5.0, 10.0);

    let min_length = if acc.min_length == usize::MAX {
        0
    } else {
        acc.min_length
    };
    let sequence_length = if min_length == acc.max_length {
        min_length.to_string()
    } else {
        format!("{min_length}-{}", acc.max_length)
    };
    let mean_length = if acc.total_sequences == 0 {
        0.0
    } else {
        acc.total_bases as f64 / acc.total_sequences as f64
    };
    let gc_percent = if acc.total_bases == 0 {
        0.0
    } else {
        acc.total_gc as f64 * 100.0 / acc.total_bases as f64
    };

    Ok(Report {
        schema_version: 1,
        program: "RustQC",
        version: env!("CARGO_PKG_VERSION"),
        generated_at_unix: SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_secs(),
        basic_statistics: Module {
            status: Status::Pass,
            data: BasicStatistics {
                filename: acc.filename,
                file_type: if acc.has_quality {
                    "Conventional base calls"
                } else {
                    "Sequence only"
                },
                encoding,
                total_sequences: acc.total_sequences,
                total_bases: acc.total_bases,
                sequences_flagged_as_poor_quality: acc.poor_quality_sequences,
                sequence_length,
                min_length,
                max_length: acc.max_length,
                mean_length,
                gc_percent,
            },
        },
        per_base_sequence_quality: Module {
            status: per_base_quality_status,
            data: per_base_quality,
        },
        per_tile_sequence_quality: Module {
            status: tile_status,
            data: tile_data,
        },
        per_sequence_quality_scores: Module {
            status: sequence_quality_status,
            data: sequence_quality,
        },
        per_base_sequence_content: Module {
            status: base_content_status,
            data: base_content,
        },
        per_sequence_gc_content: Module {
            status: gc_status,
            data: gc_content,
        },
        per_base_n_content: Module {
            status: n_status,
            data: n_content,
        },
        sequence_length_distribution: Module {
            status: length_status,
            data: lengths,
        },
        sequence_duplication_levels: Module {
            status: duplication_status,
            data: duplication,
        },
        overrepresented_sequences: Module {
            status: overrepresented_status,
            data: overrepresented,
        },
        adapter_content: Module {
            status: adapter_status,
            data: adapters,
        },
    })
}

fn quality_encoding(acc: &Accumulator, requested: PhredOffset) -> Result<(u8, String), String> {
    if !acc.has_quality {
        return Ok((33, "No quality scores".into()));
    }
    let offset = match requested {
        PhredOffset::Auto if acc.min_quality_ascii < 64 => 33,
        PhredOffset::Auto => 64,
        PhredOffset::Phred33 => 33,
        PhredOffset::Phred64 => 64,
    };
    if acc.min_quality_ascii < offset {
        return Err(format!(
            "{}: quality byte {} is below forced Phred+{} offset",
            acc.filename, acc.min_quality_ascii, offset
        ));
    }
    let name = if offset == 33 {
        "Sanger / Illumina 1.9"
    } else if acc.min_quality_ascii == 65 {
        "Illumina 1.3"
    } else {
        "Illumina 1.5"
    };
    Ok((offset, name.into()))
}

fn finish_per_base_quality(
    acc: &Accumulator,
    groups: &[PositionGroup],
    offset: u8,
) -> Vec<BaseQualityPoint> {
    groups
        .iter()
        .map(|group| {
            let positions = observed_positions(acc, group);
            let mean = average_values(
                positions
                    .clone()
                    .filter(|position| position.quality_count > 0)
                    .map(|position| {
                        position.quality_sum_ascii as f64 / position.quality_count as f64
                            - offset as f64
                    }),
            )
            .unwrap_or(0.0);
            let percentile = |percent| {
                average_values(
                    positions
                        .clone()
                        .filter(|position| position.quality_count > 100)
                        .map(|position| quality_percentile(position, offset, percent)),
                )
            };
            BaseQualityPoint {
                base: group.label(),
                position: group.midpoint(),
                mean,
                median: percentile(50),
                lower_quartile: percentile(25),
                upper_quartile: percentile(75),
                percentile_10: percentile(10),
                percentile_90: percentile(90),
            }
        })
        .collect()
}

fn quality_status(points: &[BaseQualityPoint], has_quality: bool) -> Status {
    if !has_quality {
        return Status::Skip;
    }
    let mut status = Status::Pass;
    for point in points {
        let (Some(lower), Some(median)) = (point.lower_quartile, point.median) else {
            continue;
        };
        if lower < 5.0 || median < 20.0 {
            return Status::Fail;
        }
        if lower < 10.0 || median < 25.0 {
            status = Status::Warn;
        }
    }
    status
}

fn quality_percentile(position: &PositionAccumulator, offset: u8, percentile: u64) -> f64 {
    let target = position.quality_count * percentile / 100;
    let mut cumulative = 0;
    for (index, count) in position.quality_counts.iter().enumerate() {
        cumulative += count;
        if cumulative >= target {
            return (index + QUALITY_MIN_ASCII).saturating_sub(offset as usize) as f64;
        }
    }
    0.0
}

fn finish_tile_quality(acc: &Accumulator, groups: &[PositionGroup]) -> TileQualityData {
    if !acc.tile_supported || acc.tile_counts.is_empty() {
        return TileQualityData::default();
    }
    let mut tiles: Vec<_> = acc.tile_counts.keys().copied().collect();
    tiles.sort_unstable();
    let mut max_deviation = 0.0f64;
    let deviations = tiles
        .iter()
        .map(|tile| {
            let tile_positions = &acc.tile_counts[tile];
            groups
                .iter()
                .map(|group| {
                    let mut total = 0.0;
                    let mut count = 0;
                    for raw_index in group.start - 1..group.end {
                        let Some(tile_value) = tile_positions.get(raw_index) else {
                            continue;
                        };
                        let Some(global) = acc.positions.get(raw_index) else {
                            continue;
                        };
                        if tile_value.count > 0 && global.quality_count > 0 {
                            total += tile_value.sum as f64 / tile_value.count as f64
                                - global.quality_sum_ascii as f64 / global.quality_count as f64;
                            count += 1;
                        }
                    }
                    let value = if count == 0 {
                        0.0
                    } else {
                        total / count as f64
                    };
                    max_deviation = max_deviation.max(value.abs());
                    value
                })
                .collect()
        })
        .collect();
    TileQualityData {
        tiles,
        positions: groups.iter().map(PositionGroup::label).collect(),
        position_midpoints: groups.iter().map(PositionGroup::midpoint).collect(),
        deviations,
        max_deviation,
    }
}

fn finish_sequence_quality(acc: &Accumulator, offset: u8) -> SequenceQualityData {
    let mut distribution = Vec::new();
    let mut mode = 0;
    let mut max_count = 0;
    for (ascii, count) in acc.sequence_quality_ascii.iter().copied().enumerate() {
        if count == 0 {
            continue;
        }
        let quality = ascii.saturating_sub(offset as usize) as u64;
        if count > max_count {
            max_count = count;
            mode = quality as i32;
        }
        distribution.push(CountPoint {
            value: quality,
            label: quality.to_string(),
            count: count as f64,
        });
    }
    SequenceQualityData {
        distribution,
        most_frequent_score: mode,
    }
}

fn finish_base_content(acc: &Accumulator, groups: &[PositionGroup]) -> Vec<BaseContentPoint> {
    groups
        .iter()
        .map(|group| {
            let mut counts = [0u64; 4];
            for position in observed_positions(acc, group) {
                for (total, value) in counts.iter_mut().zip(position.bases) {
                    *total += value;
                }
            }
            let total: u64 = counts.iter().sum();
            let percent = |value| {
                if total == 0 {
                    0.0
                } else {
                    value as f64 * 100.0 / total as f64
                }
            };
            BaseContentPoint {
                base: group.label(),
                position: group.midpoint(),
                a: percent(counts[0]),
                c: percent(counts[1]),
                g: percent(counts[2]),
                t: percent(counts[3]),
            }
        })
        .collect()
}

fn finish_n_content(acc: &Accumulator, groups: &[PositionGroup]) -> Vec<PositionValue> {
    groups
        .iter()
        .map(|group| {
            let mut ambiguous = 0u64;
            let mut total = 0u64;
            for position in observed_positions(acc, group) {
                ambiguous += position.bases[4];
                total += position.bases.iter().sum::<u64>();
            }
            PositionValue {
                base: group.label(),
                position: group.midpoint(),
                value: if total == 0 {
                    0.0
                } else {
                    ambiguous as f64 * 100.0 / total as f64
                },
            }
        })
        .collect()
}

fn finish_gc_content(distribution: &[f64; 101]) -> GcContentData {
    let total: f64 = distribution.iter().sum();
    if total <= 1.0 {
        return GcContentData {
            distribution: distribution
                .iter()
                .enumerate()
                .map(|(gc, count)| GcPoint {
                    gc_percent: gc as u8,
                    count: *count,
                    theoretical: *count,
                })
                .collect(),
            deviation_percent: 0.0,
        };
    }
    let mut first_mode = 0;
    let mut mode_value = 0.0;
    for (index, value) in distribution.iter().copied().enumerate() {
        if value > mode_value {
            first_mode = index;
            mode_value = value;
        }
    }
    let cutoff = distribution[first_mode] * 0.9;
    let mut mode_sum = 0.0;
    let mut mode_count = 0usize;
    let mut fell_off_top = true;
    for (index, count) in distribution.iter().enumerate().skip(first_mode) {
        if *count > cutoff {
            mode_sum += index as f64;
            mode_count += 1;
        } else {
            fell_off_top = false;
            break;
        }
    }
    let mut fell_off_bottom = true;
    for index in (0..first_mode).rev() {
        if distribution[index] > cutoff {
            mode_sum += index as f64;
            mode_count += 1;
        } else {
            fell_off_bottom = false;
            break;
        }
    }
    let mode = if fell_off_bottom || fell_off_top || mode_count == 0 {
        first_mode as f64
    } else {
        mode_sum / mode_count as f64
    };
    let variance = distribution
        .iter()
        .enumerate()
        .map(|(index, count)| (index as f64 - mode).powi(2) * count)
        .sum::<f64>()
        / (total - 1.0);
    let stdev = variance.sqrt();
    let mut deviation = 0.0;
    let points = distribution
        .iter()
        .enumerate()
        .map(|(gc, observed)| {
            let theoretical = if stdev > 0.0 {
                let density = (-(gc as f64 - mode).powi(2) / (2.0 * stdev * stdev)).exp()
                    / ((2.0 * std::f64::consts::PI).sqrt() * stdev);
                density * total
            } else if gc == first_mode {
                total
            } else {
                0.0
            };
            deviation += (theoretical - observed).abs();
            GcPoint {
                gc_percent: gc as u8,
                count: *observed,
                theoretical,
            }
        })
        .collect();
    GcContentData {
        distribution: points,
        deviation_percent: deviation * 100.0 / total,
    }
}

fn finish_lengths(lengths: &BTreeMap<usize, u64>) -> Vec<CountPoint> {
    lengths
        .iter()
        .map(|(length, count)| CountPoint {
            value: *length as u64,
            label: length.to_string(),
            count: *count as f64,
        })
        .collect()
}

fn finish_duplication(
    sequence_counts: &HashMap<Vec<u8>, u64>,
    count_at_limit: u64,
    total_count: u64,
) -> DuplicationData {
    let mut collated: HashMap<u64, u64> = HashMap::new();
    for count in sequence_counts.values() {
        *collated.entry(*count).or_default() += 1;
    }
    let mut percentages = [0.0f64; 16];
    let mut deduplicated_total = 0.0;
    let mut raw_total = 0.0;
    let mut collated: Vec<_> = collated.into_iter().collect();
    collated.sort_unstable_by_key(|(duplication_level, _)| *duplication_level);
    for (duplication_level, observations) in collated {
        let corrected = corrected_duplication_count(
            count_at_limit,
            total_count,
            duplication_level,
            observations,
        );
        deduplicated_total += corrected;
        raw_total += corrected * duplication_level as f64;
        let slot = duplication_slot(duplication_level);
        percentages[slot] += corrected * duplication_level as f64;
    }
    let labels = [
        "1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k", ">5k",
        ">10k+",
    ];
    if raw_total > 0.0 {
        for value in &mut percentages {
            *value = *value * 100.0 / raw_total;
        }
    }
    DuplicationData {
        total_deduplicated_percent: if raw_total == 0.0 {
            100.0
        } else {
            deduplicated_total * 100.0 / raw_total
        },
        levels: labels
            .iter()
            .zip(percentages)
            .map(|(label, percentage)| DuplicationPoint {
                level: (*label).to_string(),
                percentage_of_total: percentage,
            })
            .collect(),
    }
}

fn corrected_duplication_count(
    count_at_limit: u64,
    total_count: u64,
    duplication_level: u64,
    observations: u64,
) -> f64 {
    if count_at_limit == total_count
        || total_count.saturating_sub(observations) < count_at_limit
        || observations == 0
    {
        return observations as f64;
    }
    let limit = 1.0 - observations as f64 / (observations as f64 + 0.01);
    let mut probability_not_seen = 1.0;
    for index in 0..count_at_limit {
        let remaining = total_count - index;
        probability_not_seen *=
            remaining.saturating_sub(duplication_level) as f64 / remaining as f64;
        if probability_not_seen < limit {
            probability_not_seen = 0.0;
            break;
        }
    }
    observations as f64 / (1.0 - probability_not_seen)
}

fn duplication_slot(level: u64) -> usize {
    let zero_based = level.saturating_sub(1);
    match zero_based {
        0..=8 => zero_based as usize,
        9..=49 => 9,
        50..=99 => 10,
        100..=499 => 11,
        500..=999 => 12,
        1_000..=4_999 => 13,
        5_000..=9_999 => 14,
        _ => 15,
    }
}

fn finish_overrepresented(
    sequence_counts: &HashMap<Vec<u8>, u64>,
    total_count: u64,
    contaminants: &[NamedSequence],
) -> Vec<OverrepresentedSequence> {
    if total_count == 0 {
        return Vec::new();
    }
    let mut sequences: Vec<_> = sequence_counts
        .iter()
        .filter_map(|(sequence, count)| {
            let percentage = *count as f64 * 100.0 / total_count as f64;
            (percentage > OVERREP_WARN_PERCENT).then(|| OverrepresentedSequence {
                sequence: String::from_utf8_lossy(sequence).into_owned(),
                count: *count,
                percentage,
                possible_source: possible_source(sequence, contaminants),
            })
        })
        .collect();
    sequences.sort_by(|left, right| {
        right
            .count
            .cmp(&left.count)
            .then_with(|| left.sequence.cmp(&right.sequence))
    });
    sequences
}

fn finish_adapters(
    acc: &Accumulator,
    config: &AnalysisConfig,
    _sequence_groups: &[PositionGroup],
) -> AdapterContentData {
    let longest_adapter = config
        .adapters
        .iter()
        .map(|adapter| adapter.sequence.len())
        .max()
        .unwrap_or(0);
    let position_count = acc
        .max_length
        .checked_sub(longest_adapter)
        .map_or(0, |length| length + 1);
    let groups = make_base_groups(
        position_count,
        config.no_group,
        config.min_length.min(position_count),
    );
    let mut max_content = 0.0f64;
    let series = config
        .adapters
        .iter()
        .enumerate()
        .map(|(adapter_index, adapter)| {
            let mut cumulative = vec![0.0; position_count];
            let mut running = 0u64;
            for (position, value) in cumulative.iter_mut().enumerate() {
                running += acc.adapter_starts[adapter_index]
                    .get(position)
                    .copied()
                    .unwrap_or(0);
                *value = if acc.analyzed_sequences == 0 {
                    0.0
                } else {
                    running as f64 * 100.0 / acc.analyzed_sequences as f64
                };
            }
            let percentages = groups
                .iter()
                .map(|group| {
                    let values = &cumulative[group.start - 1..group.end];
                    let value = values.iter().sum::<f64>() / values.len() as f64;
                    max_content = max_content.max(value);
                    value
                })
                .collect();
            AdapterSeries {
                name: adapter.name.clone(),
                percentages,
            }
        })
        .collect();
    AdapterContentData {
        positions: groups.iter().map(PositionGroup::label).collect(),
        position_midpoints: groups.iter().map(PositionGroup::midpoint).collect(),
        series,
        max_content,
    }
}

fn observed_positions<'a>(
    acc: &'a Accumulator,
    group: &PositionGroup,
) -> impl Iterator<Item = &'a PositionAccumulator> + Clone {
    let start = group.start.saturating_sub(1).min(acc.positions.len());
    let end = group.end.min(acc.positions.len());
    acc.positions[start..end].iter()
}

fn average_values(values: impl Iterator<Item = f64>) -> Option<f64> {
    let (sum, count) = values.fold((0.0, 0usize), |(sum, count), value| {
        (sum + value, count + 1)
    });
    (count > 0).then_some(sum / count as f64)
}

fn gc_claim_counts(length: usize) -> [u32; 101] {
    let mut claims = [0u32; 101];
    for gc_count in 0..=length {
        let (low, high) = gc_percentage_range(gc_count, length);
        for claim_count in claims.iter_mut().take(high + 1).skip(low) {
            *claim_count += 1;
        }
    }
    claims
}

fn gc_percentage_range(gc_count: usize, length: usize) -> (usize, usize) {
    let low_count = (gc_count as f64 - 0.5).max(0.0);
    let high_count = (gc_count as f64 + 0.5).min(length as f64);
    let low = (low_count * 100.0 / length as f64).round() as usize;
    let high = (high_count * 100.0 / length as f64).round() as usize;
    (low.min(100), high.min(100))
}

fn find_subslice_case_insensitive(sequence: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || sequence.len() < needle.len() {
        return None;
    }
    sequence.windows(needle.len()).position(|window| {
        window
            .iter()
            .zip(needle)
            .all(|(left, right)| left.to_ascii_uppercase() == *right)
    })
}

fn is_casava_filtered(id: &[u8]) -> bool {
    id.split(|byte| byte.is_ascii_whitespace())
        .nth(1)
        .and_then(|suffix| suffix.split(|byte| *byte == b':').nth(1))
        .is_some_and(|filter| filter == b"Y")
}

fn parse_tile(id: &[u8]) -> Option<u32> {
    let first_field = id
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or(id);
    let fields: Vec<_> = first_field.split(|byte| *byte == b':').collect();
    let tile = if fields.len() >= 7 {
        fields.get(4)?
    } else if fields.len() >= 5 {
        fields.get(2)?
    } else {
        return None;
    };
    std::str::from_utf8(tile).ok()?.parse().ok()
}

fn input_filename(input: &str) -> String {
    if input == "-" {
        "stdin".into()
    } else if let Some(name) = input.strip_prefix("stdin:") {
        if name.is_empty() {
            "stdin".into()
        } else {
            name.to_string()
        }
    } else {
        Path::new(input)
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| input.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_modern_and_legacy_tile_ids() {
        assert_eq!(
            parse_tile(b"K00271:89:HHWWNBBXX:2:1101:9749:1086 1:N:0:A"),
            Some(1101)
        );
        assert_eq!(parse_tile(b"HWUSI-EAS493_0001:2:7:1000:16900#0/1"), Some(7));
        assert_eq!(parse_tile(b"plain-read"), None);
    }

    #[test]
    fn detects_casava_filter_flag() {
        assert!(is_casava_filtered(b"id 1:Y:0:INDEX"));
        assert!(!is_casava_filtered(b"id 1:N:0:INDEX"));
    }

    #[test]
    fn gc_model_balances_percentage_bins() {
        let claims = gc_claim_counts(50);
        let total: f64 = (0..=50)
            .map(|gc_count| {
                let (low, high) = gc_percentage_range(gc_count, 50);
                (low..=high)
                    .map(|percent| 1.0 / claims[percent] as f64)
                    .sum::<f64>()
            })
            .sum();
        assert!((total - 101.0).abs() < 1e-9);
    }

    #[test]
    fn parallel_batches_preserve_order_sensitive_statistics() {
        let path = std::env::temp_dir().join(format!(
            "rustqc-parallel-{}-{}.fastq",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("system clock should be valid")
                .as_nanos()
        ));
        let mut fastq = String::with_capacity((RECORDS_PER_BATCH + 19) * 100);
        for index in 0..RECORDS_PER_BATCH + 19 {
            let sequence = if index % 3 == 0 {
                "ACGTACGT"
            } else {
                "GGGGTTTT"
            };
            fastq.push_str(&format!(
                "@K00271:89:HHWWNBBXX:2:1101:9749:{index} 1:N:0:A\n{sequence}\n+\nIIIIIIII\n"
            ));
        }
        std::fs::write(&path, fastq).expect("fixture should be writable");
        let config = AnalysisConfig {
            no_group: false,
            nofilter: false,
            dup_length: 50,
            min_length: 0,
            phred_offset: PhredOffset::Phred33,
            adapters: crate::metrics::default_adapters(),
            contaminants: crate::metrics::default_contaminants(),
        };
        let sequential = analyze_input(path.to_str().expect("utf-8 path"), &config)
            .expect("sequential analysis should succeed");
        let parallel =
            analyze_input_with_workers(path.to_str().expect("utf-8 path"), &config, 4, false)
                .expect("parallel analysis should succeed");
        assert_eq!(
            sequential.basic_statistics.data.total_sequences,
            parallel.basic_statistics.data.total_sequences
        );
        assert_eq!(
            sequential.sequence_duplication_levels.status,
            parallel.sequence_duplication_levels.status
        );
        assert_eq!(
            sequential.per_tile_sequence_quality.data.tiles,
            parallel.per_tile_sequence_quality.data.tiles
        );
        assert_eq!(
            sequential.adapter_content.data.positions,
            parallel.adapter_content.data.positions
        );
        let _ = std::fs::remove_file(path);
    }
}
