use crate::cli::{Args, OutputFormat};
use crate::html;
use crate::report::{Report, Status};
use std::fmt::Write as FmtWrite;
use std::fs::{self, File};
use std::io::Write as IoWrite;
use std::path::{Path, PathBuf};
use zip::write::SimpleFileOptions;
use zip::{CompressionMethod, ZipWriter};

pub fn write_artifacts(report: &Report, args: &Args) -> Result<Vec<PathBuf>, String> {
    fs::create_dir_all(&args.outdir)
        .map_err(|error| format!("failed to create {}: {error}", args.outdir.display()))?;
    if let Some(path) = &args.output {
        let contents = match args.format {
            OutputFormat::Html => html::render_html(report)?,
            OutputFormat::Json => serde_json::to_string_pretty(report)
                .map_err(|error| format!("failed to serialize JSON: {error}"))?,
            OutputFormat::Text => render_fastqc_data(report),
            OutputFormat::All => unreachable!("validated by Args::all_inputs"),
        };
        write_file(path, contents.as_bytes())?;
        return Ok(vec![path.clone()]);
    }

    let stem = report_stem(&report.basic_statistics.data.filename);
    let mut paths = Vec::new();
    let html_report = matches!(args.format, OutputFormat::Html | OutputFormat::All)
        .then(|| html::render_html(report))
        .transpose()?;
    let summary = render_summary(report);
    let data = render_fastqc_data(report);

    if let Some(html) = &html_report {
        let path = args.outdir.join(format!("{stem}_rustqc.html"));
        write_file(&path, html.as_bytes())?;
        paths.push(path);
        if !args.nozip {
            let path = args.outdir.join(format!("{stem}_rustqc.zip"));
            write_zip(&path, &stem, html, &summary, &data)?;
            paths.push(path);
        }
    }

    if matches!(args.format, OutputFormat::Json | OutputFormat::All) {
        let path = args.outdir.join(format!("{stem}_rustqc.json"));
        let json = serde_json::to_vec_pretty(report)
            .map_err(|error| format!("failed to serialize JSON: {error}"))?;
        write_file(&path, &json)?;
        paths.push(path);
    }

    if matches!(args.format, OutputFormat::Text | OutputFormat::All) {
        let data_path = args.outdir.join(format!("{stem}_fastqc_data.txt"));
        write_file(&data_path, data.as_bytes())?;
        paths.push(data_path);
        let summary_path = args.outdir.join(format!("{stem}_summary.txt"));
        write_file(&summary_path, summary.as_bytes())?;
        paths.push(summary_path);
    }

    if args.extract {
        let directory = args.outdir.join(format!("{stem}_rustqc"));
        fs::create_dir_all(&directory)
            .map_err(|error| format!("failed to create {}: {error}", directory.display()))?;
        write_file(&directory.join("summary.txt"), summary.as_bytes())?;
        write_file(&directory.join("fastqc_data.txt"), data.as_bytes())?;
        if let Some(html) = &html_report {
            write_file(&directory.join("fastqc_report.html"), html.as_bytes())?;
        }
        paths.push(directory);
    }
    Ok(paths)
}

pub fn render_summary(report: &Report) -> String {
    let filename = sanitize_field(&report.basic_statistics.data.filename);
    let mut summary = String::with_capacity(768);
    for (name, status) in report.summaries() {
        if status != Status::Skip {
            let _ = writeln!(summary, "{status}\t{name}\t{filename}");
        }
    }
    summary
}

pub fn render_fastqc_data(report: &Report) -> String {
    let mut output = String::with_capacity(64 * 1024);
    let basic = &report.basic_statistics.data;
    let _ = writeln!(output, "##FastQC\t0.12.1");
    let _ = writeln!(output, "##RustQC\t{}", report.version);
    begin_module(
        &mut output,
        "Basic Statistics",
        report.basic_statistics.status,
    );
    output.push_str("#Measure\tValue\n");
    let _ = writeln!(output, "Filename\t{}", sanitize_field(&basic.filename));
    let _ = writeln!(output, "File type\t{}", basic.file_type);
    let _ = writeln!(output, "Encoding\t{}", basic.encoding);
    let _ = writeln!(output, "Total Sequences\t{}", basic.total_sequences);
    let _ = writeln!(
        output,
        "Total Bases\t{}",
        format_base_count(basic.total_bases)
    );
    let _ = writeln!(
        output,
        "Sequences flagged as poor quality\t{}",
        basic.sequences_flagged_as_poor_quality
    );
    let _ = writeln!(output, "Sequence length\t{}", basic.sequence_length);
    let _ = writeln!(output, "%GC\t{:.0}", basic.gc_percent);
    end_module(&mut output);

    if report.per_base_sequence_quality.status != Status::Skip {
        begin_module(
            &mut output,
            "Per base sequence quality",
            report.per_base_sequence_quality.status,
        );
        output.push_str(
            "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n",
        );
        for point in &report.per_base_sequence_quality.data {
            let _ = writeln!(
                output,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                point.base,
                point.mean,
                display_option(point.median),
                display_option(point.lower_quartile),
                display_option(point.upper_quartile),
                display_option(point.percentile_10),
                display_option(point.percentile_90)
            );
        }
        end_module(&mut output);
    }

    let tile = &report.per_tile_sequence_quality;
    if tile.status != Status::Skip {
        begin_module(&mut output, "Per tile sequence quality", tile.status);
        output.push_str("#Tile\tBase\tMean\n");
        for (tile_id, deviations) in tile.data.tiles.iter().zip(&tile.data.deviations) {
            for (base, deviation) in tile.data.positions.iter().zip(deviations) {
                let _ = writeln!(output, "{tile_id}\t{base}\t{deviation}");
            }
        }
        end_module(&mut output);
    }

    let quality = &report.per_sequence_quality_scores;
    if quality.status != Status::Skip {
        begin_module(&mut output, "Per sequence quality scores", quality.status);
        output.push_str("#Quality\tCount\n");
        for point in &quality.data.distribution {
            let _ = writeln!(output, "{}\t{}", point.label, point.count);
        }
        end_module(&mut output);
    }

    begin_module(
        &mut output,
        "Per base sequence content",
        report.per_base_sequence_content.status,
    );
    output.push_str("#Base\tG\tA\tT\tC\n");
    for point in &report.per_base_sequence_content.data {
        let _ = writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}",
            point.base, point.g, point.a, point.t, point.c
        );
    }
    end_module(&mut output);

    begin_module(
        &mut output,
        "Per sequence GC content",
        report.per_sequence_gc_content.status,
    );
    output.push_str("#GC Content\tCount\n");
    for point in &report.per_sequence_gc_content.data.distribution {
        let _ = writeln!(output, "{}\t{}", point.gc_percent, point.count);
    }
    end_module(&mut output);

    begin_module(
        &mut output,
        "Per base N content",
        report.per_base_n_content.status,
    );
    output.push_str("#Base\tN-Count\n");
    for point in &report.per_base_n_content.data {
        let _ = writeln!(output, "{}\t{}", point.base, point.value);
    }
    end_module(&mut output);

    begin_module(
        &mut output,
        "Sequence Length Distribution",
        report.sequence_length_distribution.status,
    );
    output.push_str("#Length\tCount\n");
    for point in &report.sequence_length_distribution.data {
        let _ = writeln!(output, "{}\t{}", point.label, point.count);
    }
    end_module(&mut output);

    let duplication = &report.sequence_duplication_levels;
    begin_module(
        &mut output,
        "Sequence Duplication Levels",
        duplication.status,
    );
    let _ = writeln!(
        output,
        "#Total Deduplicated Percentage\t{}",
        duplication.data.total_deduplicated_percent
    );
    output.push_str("#Duplication Level\tPercentage of total\n");
    for point in &duplication.data.levels {
        let _ = writeln!(output, "{}\t{}", point.level, point.percentage_of_total);
    }
    end_module(&mut output);

    let overrepresented = &report.overrepresented_sequences;
    begin_module(
        &mut output,
        "Overrepresented sequences",
        overrepresented.status,
    );
    output.push_str("#Sequence\tCount\tPercentage\tPossible Source\n");
    for sequence in &overrepresented.data {
        let _ = writeln!(
            output,
            "{}\t{}\t{}\t{}",
            sanitize_field(&sequence.sequence),
            sequence.count,
            sequence.percentage,
            sanitize_field(&sequence.possible_source)
        );
    }
    end_module(&mut output);

    let adapters = &report.adapter_content;
    begin_module(&mut output, "Adapter Content", adapters.status);
    output.push_str("#Position");
    for series in &adapters.data.series {
        let _ = write!(output, "\t{}", sanitize_field(&series.name));
    }
    output.push('\n');
    for (index, position) in adapters.data.positions.iter().enumerate() {
        output.push_str(position);
        for series in &adapters.data.series {
            let _ = write!(output, "\t{}", series.percentages[index]);
        }
        output.push('\n');
    }
    end_module(&mut output);
    output
}

fn write_zip(path: &Path, stem: &str, html: &str, summary: &str, data: &str) -> Result<(), String> {
    let file = File::create(path)
        .map_err(|error| format!("failed to create {}: {error}", path.display()))?;
    let mut archive = ZipWriter::new(file);
    let options = SimpleFileOptions::default()
        .compression_method(CompressionMethod::Deflated)
        .compression_level(Some(1))
        .unix_permissions(0o644);
    let root = format!("{stem}_fastqc");
    for (name, contents) in [
        ("summary.txt", summary),
        ("fastqc_data.txt", data),
        ("fastqc_report.html", html),
    ] {
        archive
            .start_file(format!("{root}/{name}"), options)
            .map_err(|error| format!("failed to write {}: {error}", path.display()))?;
        archive
            .write_all(contents.as_bytes())
            .map_err(|error| format!("failed to write {}: {error}", path.display()))?;
    }
    archive
        .finish()
        .map_err(|error| format!("failed to finish {}: {error}", path.display()))?;
    Ok(())
}

fn write_file(path: &Path, contents: &[u8]) -> Result<(), String> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)
            .map_err(|error| format!("failed to create {}: {error}", parent.display()))?;
    }
    fs::write(path, contents)
        .map_err(|error| format!("failed to write {}: {error}", path.display()))
}

fn begin_module(output: &mut String, name: &str, status: Status) {
    let _ = writeln!(
        output,
        ">>{name}\t{}",
        status.to_string().to_ascii_lowercase()
    );
}

fn end_module(output: &mut String) {
    output.push_str(">>END_MODULE\n");
}

fn display_option(value: Option<f64>) -> String {
    value.map_or_else(|| "NaN".into(), |value| value.to_string())
}

fn sanitize_field(value: &str) -> String {
    value.replace(['\t', '\r', '\n'], " ")
}

fn format_base_count(count: u64) -> String {
    if count >= 1_000_000_000 {
        format!("{:.0} Gbp", count as f64 / 1_000_000_000.0)
    } else if count >= 1_000_000 {
        format!("{:.0} Mbp", count as f64 / 1_000_000.0)
    } else if count >= 1_000 {
        format!("{:.0} kbp", count as f64 / 1_000.0)
    } else {
        count.to_string()
    }
}

fn report_stem(filename: &str) -> String {
    let lower = filename.to_ascii_lowercase();
    for suffix in [
        ".fastq.gz",
        ".fq.gz",
        ".fasta.gz",
        ".fa.gz",
        ".fastq.bz2",
        ".fq.bz2",
        ".fastq",
        ".fq",
        ".fasta",
        ".fa",
        ".txt",
    ] {
        if lower.ends_with(suffix) {
            return filename[..filename.len() - suffix.len()].to_string();
        }
    }
    Path::new(filename)
        .file_stem()
        .map(|stem| stem.to_string_lossy().into_owned())
        .filter(|stem| !stem.is_empty())
        .unwrap_or_else(|| "rustqc".into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_compound_fastq_suffixes() {
        assert_eq!(report_stem("sample.fastq.gz"), "sample");
        assert_eq!(report_stem("sample.fq"), "sample");
    }

    #[test]
    fn formats_fastqc_base_units() {
        assert_eq!(format_base_count(19_059_150), "19 Mbp");
    }
}
