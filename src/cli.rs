use clap::Parser;

/// RustQC: A lightweight, fast quality control tool for FASTQ files.
#[derive(Parser, Debug)]
#[command(author = "Shakil Rafi", version, about = "FastQC in Rust", long_about = None)]
pub struct Args {
    /// Input FASTQ file (can be .gz)
    #[arg(short, long)]
    pub input: String,

    /// Output format: json, markdown, etc.
    #[arg(short, long, default_value = "json")]
    pub format: String,

    /// Optional output file path
    #[arg(short, long)]
    pub output: Option<String>,
}
