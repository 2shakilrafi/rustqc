use needletail::parse_fastx_file;
use crate::report::Report;
use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use indicatif::{ProgressBar, ProgressStyle};

pub fn analyze_fastq(path: &str) -> Report {
    let file_size = std::fs::metadata(path).expect("Failed to get metadata").len();

    let progress = ProgressBar::new(file_size);
    progress.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({eta})")
        .unwrap()
        .progress_chars("█▓░"));

    let mut reader = parse_fastx_file(path).expect("Failed to open FASTQ file");

    let mut total_reads = 0;
    let mut total_bases = 0;
    let mut total_gc = 0;

    let mut pos_quality: Vec<Vec<u8>> = Vec::new();
    let mut length_counts: HashMap<usize, usize> = HashMap::new();
    let mut read_length_histogram = HashMap::with_capacity(1024);
    let mut gc_percent_histogram = HashMap::with_capacity(101);

    let mut base_counts_a: Vec<usize> = Vec::new();
    let mut base_counts_t: Vec<usize> = Vec::new();
    let mut base_counts_g: Vec<usize> = Vec::new();
    let mut base_counts_c: Vec<usize> = Vec::new();
    let mut base_counts_total: Vec<usize> = Vec::new();

    let mut min_length = usize::MAX;
    let mut max_length = 0;
    let mut per_sequence_quality_histogram = HashMap::with_capacity(64);

    let mut n_content_counts: Vec<usize> = Vec::new();

    let mut kmer_counts: HashMap<Vec<u8>, usize> = HashMap::new();
    let kmer_len = 10;



    while let Some(record) = reader.next() {
        let seqrec = record.expect("Invalid FASTQ record");
        let seq = seqrec.seq();
        let qual = seqrec.qual().unwrap_or(&[]);

        let len = seq.len();
        if len == 0 {
            continue;
        }

        total_reads += 1;
        total_bases += len;

        let read_gc = seq.iter().filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c').count();
        total_gc += read_gc;

        let gc_percent = ((read_gc as f64 / len as f64) * 100.0).round() as u8;
        *gc_percent_histogram.entry(gc_percent).or_insert(0) += 1;

        *length_counts.entry(len).or_insert(0) += 1;
        *read_length_histogram.entry(len).or_insert(0) += 1;

        min_length = min_length.min(len);
        max_length = max_length.max(len);

        for (i, &q) in qual.iter().enumerate() {
            if pos_quality.len() <= i {
                pos_quality.push(Vec::new());
            }
            pos_quality[i].push(q);
        }

        for (i, &b) in seq.iter().enumerate() {
            if base_counts_total.len() <= i {
                base_counts_total.resize(i + 1, 0);
                base_counts_a.resize(i + 1, 0);
                base_counts_t.resize(i + 1, 0);
                base_counts_g.resize(i + 1, 0);
                base_counts_c.resize(i + 1, 0);
                n_content_counts.resize(i + 1, 0);  // ← NEW
            }

            base_counts_total[i] += 1;

            match b.to_ascii_uppercase() {
                b'A' => base_counts_a[i] += 1,
                b'T' => base_counts_t[i] += 1,
                b'G' => base_counts_g[i] += 1,
                b'C' => base_counts_c[i] += 1,
                b'N' => n_content_counts[i] += 1,  // ← NEW
                _ => {}
            }
        }

        let n_content: Vec<f64> = n_content_counts.iter()
            .zip(&base_counts_total)
            .map(|(&n, &t)| (n as f64 / t as f64) * 100.0)
            .collect();



        if !qual.is_empty() {
            let avg_q = (qual.iter().map(|&q| q as usize).sum::<usize>() as f64 / qual.len() as f64).round() as u8;
            *per_sequence_quality_histogram.entry(avg_q).or_insert(0) += 1;
        }

        // Best effort progress estimation
        let estimated_bytes = (total_bases as f64 * 1.5) as u64;
        progress.set_position(estimated_bytes.min(file_size));

        if seq.len() >= kmer_len {
            let kmer = seq[..kmer_len].to_vec();
            *kmer_counts.entry(kmer).or_insert(0) += 1;
        }

    }

    progress.finish_with_message("Done");

    let mut overrepresented_sequences: Vec<(String, f64, usize)> = kmer_counts
    .into_iter()
    .map(|(kmer, count)| {
        let kmer_str = String::from_utf8_lossy(&kmer).into_owned();
        let percentage = (count as f64 / total_reads as f64) * 100.0;
        (kmer_str, percentage, count)
    })
    .filter(|(_, pct, _)| *pct > 0.1) // Only keep kmers that show up in >0.1% of reads
    .collect();

    overrepresented_sequences.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());


    let avg_length = total_bases as f64 / total_reads as f64;
    let mode_length = length_counts
        .into_iter()
        .max_by_key(|&(_, count)| count)
        .map(|(length, _)| length)
        .unwrap_or(0);

    let qualities_per_position: Vec<f64> = pos_quality
        .into_iter()
        .map(|vals| {
            let sum: usize = vals.iter().map(|&v| v as usize).sum();
            sum as f64 / vals.len() as f64
        })
        .collect();

    let gc_content = total_gc as f64 / total_bases as f64 * 100.0;

    let percent_a = base_counts_a.iter().zip(&base_counts_total).map(|(&a, &t)| (a as f64 / t as f64) * 100.0).collect();
    let percent_t = base_counts_t.iter().zip(&base_counts_total).map(|(&a, &t)| (a as f64 / t as f64) * 100.0).collect();
    let percent_g = base_counts_g.iter().zip(&base_counts_total).map(|(&a, &t)| (a as f64 / t as f64) * 100.0).collect();
    let percent_c = base_counts_c.iter().zip(&base_counts_total).map(|(&a, &t)| (a as f64 / t as f64) * 100.0).collect();

    let n_content: Vec<f64> = n_content_counts.iter()
            .zip(&base_counts_total)
            .map(|(&n, &t)| (n as f64 / t as f64) * 100.0)
            .collect();
    
    Report {
        filename: Path::new(path)
            .file_name()
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "unknown".into()),
        total_reads,
        avg_length,
        min_length,
        max_length,
        mode_length,
        gc_content,
        qualities_per_position,
        read_length_histogram,
        gc_percent_histogram,
        percent_a,
        percent_t,
        percent_g,
        percent_c,
        per_sequence_quality_histogram,
        n_content,
        overrepresented_sequences,
    }
}
