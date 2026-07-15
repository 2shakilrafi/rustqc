# RustQC

RustQC 0.3.6 is a fast, memory-bounded command-line quality-control analyzer for FASTQ and FASTA data. It implements the eleven active FastQC 0.12.1 analysis modules and emits self-contained reports plus FastQC-compatible data artifacts.

## Analysis modules

- Basic Statistics
- Per base sequence quality, including mean, median, quartiles, and 10th/90th percentiles
- Per tile sequence quality for modern and legacy Illumina read identifiers
- Per sequence quality scores
- Per base sequence content
- Per sequence GC content with FastQC-compatible GC interpolation and normality scoring
- Per base N content
- Sequence Length Distribution
- Sequence Duplication Levels with bounded extrapolation after 100,000 unique sequences
- Overrepresented sequences with contaminant matching
- Adapter Content for Illumina, small RNA, Nextera, PolyA, and PolyG adapters

The default warn/fail thresholds and position grouping match FastQC 0.12.1. FASTA inputs are accepted, with quality-only modules marked as skipped.

## Performance design

- One streaming pass through each input; only a small, bounded number of record batches are in flight.
- Fixed quality histograms replace per-base quality vectors.
- Duplicate tracking is bounded at 100,000 unique sequence prefixes.
- Per-tile analysis uses FastQC's first-10,000-plus-10-percent sampling rule.
- A single input uses a bounded parser/decompressor, worker, and ordered-reduction pipeline; CPU-heavy base, quality, GC, and adapter statistics run concurrently.
- Rayon schedules independent input files concurrently; the `--threads` budget is divided across inputs without oversubscribing the host.
- Interactive runs show a byte-progress bar with a terminal spinner. Use `--quiet` or `--no-progress` for log-friendly output.
- Reports use embedded uPlot assets and direct serialization, with no runtime template engine or network requests.
- Release builds use thin LTO, one codegen unit, stripped symbols, and abort-on-panic.

Memory grows with the longest observed read, capped duplication state, and sampled tile state rather than the total number of reads.

## Build

```bash
cargo build --release
```

The optimized binary is `target/release/rustqc`.

## Use

Generate a self-contained HTML report and a FastQC-compatible zip archive:

```bash
rustqc sample.fastq.gz
```

Analyze a large single file with eight analysis workers:

```bash
rustqc --threads 8 --outdir qc-results sample.fastq.gz
```

Analyze multiple files concurrently:

```bash
rustqc --threads 4 --outdir qc-results reads/*.fastq.gz
```

Generate every artifact:

```bash
rustqc --format all --extract --outdir qc-results sample.fastq.gz
```

Read uncompressed FASTQ from stdin and choose the report name:

```bash
gzip -dc sample.fastq.gz | rustqc stdin:sample
```

Use custom adapter or contaminant definitions:

```bash
rustqc --adapters adapters.tsv --contaminants contaminants.tsv sample.fastq.gz
```

Each definition file is tab-separated:

```text
Display name<TAB>ACGTACGTACGT
```

Run `rustqc --help` for the complete CLI.

## Outputs

The default run creates:

- `<sample>_rustqc.html`: self-contained interactive report
- `<sample>_rustqc.zip`: archive containing `summary.txt`, `fastqc_data.txt`, and `fastqc_report.html`

`--format json` writes a stable structured report. `--format text` writes standalone FastQC-compatible summary and data files. `--format all` writes all formats. `--extract` also writes the archive contents as a directory.

The compatibility zip uses the conventional `<sample>_fastqc/fastqc_data.txt` layout so existing pipeline collectors can discover it.

## Validation

```bash
cargo fmt --all -- --check
cargo test
cargo clippy --all-targets -- -D warnings
```

For a local performance sample with macOS peak-memory reporting:

```bash
/usr/bin/time -l target/release/rustqc --quiet --format json --output report.json sample.fastq.gz
```

RustQC is CLI-first and does not reproduce FastQC's Swing desktop interface. Its parity target is FastQC 0.12.1's active sequence-QC modules and pipeline report artifacts.
