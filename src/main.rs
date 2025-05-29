mod cli;
mod reader;
mod report;
mod html;


use clap::Parser;
use cli::Args;

fn main() {
    let args = Args::parse();
    let report = reader::analyze_fastq(&args.input);

    match args.format.as_str() {
        "html" => html::render_html(&report, args.output.as_deref()),
        _ => {
            let json = serde_json::to_string_pretty(&report).unwrap();
            if let Some(path) = args.output {
                std::fs::write(path, json).expect("Failed to write output");
            } else {
                println!("{}", json);
            }
        }
    }
}

