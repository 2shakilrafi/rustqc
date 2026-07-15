use clap::{Parser, ValueEnum};
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum OutputFormat {
    Html,
    Json,
    Text,
    All,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub enum PhredOffset {
    Auto,
    #[value(name = "33")]
    Phred33,
    #[value(name = "64")]
    Phred64,
}

/// Fast, memory-bounded quality control for FASTQ data.
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// FASTQ/FASTA inputs. Compression is detected automatically. Use '-' or stdin:name for stdin.
    #[arg(value_name = "INPUT")]
    pub inputs: Vec<String>,

    /// Legacy single-input spelling retained for rustqc 0.1 compatibility.
    #[arg(short = 'i', long = "input", value_name = "INPUT", hide = true)]
    pub legacy_input: Option<String>,

    /// Directory for generated reports.
    #[arg(short = 'o', long = "outdir", default_value = ".")]
    pub outdir: PathBuf,

    /// Exact output path for a single input and a single non-HTML format.
    #[arg(long, value_name = "PATH")]
    pub output: Option<PathBuf>,

    /// Report artifacts to create.
    #[arg(short = 'f', long, value_enum, default_value_t = OutputFormat::Html)]
    pub format: OutputFormat,

    /// Total analysis workers; a single input is processed in parallel. 0 uses available CPUs.
    #[arg(short = 't', long, default_value_t = 0)]
    pub threads: usize,

    /// Also write an extracted report directory beside the archive.
    #[arg(long)]
    pub extract: bool,

    /// Do not create the compatibility zip archive for HTML output.
    #[arg(long)]
    pub nozip: bool,

    /// Disable FastQC-style position grouping.
    #[arg(long = "nogroup", alias = "no-group")]
    pub no_group: bool,

    /// Include reads marked as failed by the Illumina CASAVA filter.
    #[arg(long)]
    pub nofilter: bool,

    /// Truncate sequences to this many bases for duplication analysis (default: 50).
    #[arg(long, default_value_t = 50)]
    pub dup_length: usize,

    /// Force quality encoding or detect it from the lowest observed quality byte.
    #[arg(long, value_enum, default_value_t = PhredOffset::Auto)]
    pub phred_offset: PhredOffset,

    /// Custom tab-separated adapter list: name<TAB>sequence.
    #[arg(long, value_name = "FILE")]
    pub adapters: Option<PathBuf>,

    /// Custom tab-separated contaminant list: name<TAB>sequence.
    #[arg(long, value_name = "FILE")]
    pub contaminants: Option<PathBuf>,

    /// Extend plots to at least this sequence length.
    #[arg(long, default_value_t = 0)]
    pub min_length: usize,

    /// Suppress progress messages.
    #[arg(short = 'q', long)]
    pub quiet: bool,

    /// Disable the interactive progress indicator.
    #[arg(long)]
    pub no_progress: bool,
}

impl Args {
    pub fn all_inputs(&self) -> Result<Vec<String>, String> {
        let mut inputs = self.inputs.clone();
        if let Some(input) = &self.legacy_input {
            inputs.push(input.clone());
        }
        if inputs.is_empty() {
            return Err("at least one INPUT is required".to_string());
        }
        let stdin_count = inputs
            .iter()
            .filter(|input| input.as_str() == "-" || input.starts_with("stdin:"))
            .count();
        if stdin_count > 1 || (stdin_count == 1 && inputs.len() != 1) {
            return Err("stdin can only be analyzed as the sole input".to_string());
        }
        if self.output.is_some() && (inputs.len() != 1 || self.format == OutputFormat::All) {
            return Err(
                "--output requires one input and cannot be combined with --format all".into(),
            );
        }
        if self.dup_length == 0 {
            return Err("--dup-length must be greater than zero".into());
        }
        Ok(inputs)
    }
}
