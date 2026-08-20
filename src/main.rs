mod cli;
mod html;
mod metrics;
mod output;
mod reader;
mod report;

use clap::Parser;
use cli::Args;
use metrics::{default_adapters, default_contaminants, load_named_sequences};
use rayon::prelude::*;
use reader::{AnalysisConfig, analyze_input_with_workers};

fn main() {
    if let Err(error) = run() {
        eprintln!("rustqc: {error}");
        std::process::exit(2);
    }
}

fn run() -> Result<(), String> {
    let args = Args::parse();
    let inputs = args.all_inputs()?;
    let adapters = args
        .adapters
        .as_deref()
        .map(load_named_sequences)
        .transpose()?
        .unwrap_or_else(default_adapters);
    let contaminants = args
        .contaminants
        .as_deref()
        .map(load_named_sequences)
        .transpose()?
        .unwrap_or_else(default_contaminants);
    let config = AnalysisConfig::new(
        args.no_group,
        args.nofilter,
        args.dup_length,
        args.min_length,
        args.phred_offset,
        adapters,
        contaminants,
    )?;
    let threads = if args.threads == 0 {
        std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1)
    } else {
        args.threads
    };
    let input_workers = threads.min(inputs.len().max(1));
    // Each Rayon input task owns the parser and ordered reducer; reserve that core
    // before allocating CPU-only batch workers so concurrent inputs stay in budget.
    let workers_per_input = ((threads.saturating_sub(input_workers)) / input_workers).max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(input_workers)
        .thread_name(|index| format!("rustqc-{index}"))
        .build()
        .map_err(|error| format!("failed to initialize worker pool: {error}"))?;

    if !args.quiet {
        eprintln!(
            "RustQC {}: analyzing {} input(s) with {} worker(s)",
            env!("CARGO_PKG_VERSION"),
            inputs.len(),
            threads
        );
    }
    let results = pool.install(|| {
        inputs
            .par_iter()
            .map(|input| {
                (
                    input,
                    analyze_input_with_workers(
                        input,
                        &config,
                        workers_per_input,
                        !args.quiet && !args.no_progress && inputs.len() == 1,
                    ),
                )
            })
            .collect::<Vec<_>>()
    });

    let mut failures = Vec::new();
    for (input, result) in results {
        match result {
            Ok(report) => match output::write_artifacts(&report, &args) {
                Ok(paths) => {
                    if !args.quiet {
                        eprintln!(
                            "Completed {input}: {}",
                            paths
                                .iter()
                                .map(|path| path.display().to_string())
                                .collect::<Vec<_>>()
                                .join(", ")
                        );
                    }
                }
                Err(error) => failures.push(format!("{input}: {error}")),
            },
            Err(error) => failures.push(error),
        }
    }
    if failures.is_empty() {
        Ok(())
    } else {
        Err(failures.join("\n"))
    }
}
