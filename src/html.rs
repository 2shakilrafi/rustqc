use tera::{Tera, Context};
use crate::report::Report;
use std::fs::File;
use std::io::Write;

pub fn render_html(report: &Report, out_path: Option<&str>) {
    let tera = match Tera::new("templates/*.html.tera") {
        Ok(t) => t,
        Err(e) => {
            eprintln!("Template parsing error: {}", e);
            return;
        }
    };

    let mut context = Context::new();
    context.insert("filename", &report.filename);
    context.insert("total_reads", &report.total_reads);
    context.insert("avg_length", &report.avg_length);
    context.insert("min_length", &report.min_length);
    context.insert("max_length", &report.max_length);
    context.insert("mode_length", &report.mode_length);
    context.insert("gc_content", &report.gc_content);
    context.insert("qualities_per_position", &report.qualities_per_position);
    context.insert("read_length_histogram", &report.read_length_histogram);
    context.insert("gc_percent_histogram", &report.gc_percent_histogram);

    // ✅ New base composition vectors
    context.insert("percent_a", &report.percent_a);
    context.insert("percent_t", &report.percent_t);
    context.insert("percent_g", &report.percent_g);
    context.insert("percent_c", &report.percent_c);

    context.insert("per_sequence_quality_histogram", &report.per_sequence_quality_histogram);

    context.insert("n_content", &report.n_content);

    let now = chrono::Local::now();
    let username = whoami::username();

    context.insert("datetime", &now.format("%Y-%m-%d %H:%M:%S").to_string());
    context.insert("username", &username);
    context.insert("overrepresented_sequences", &report.overrepresented_sequences);




    let rendered = match tera.render("report.html.tera", &context) {
        Ok(html) => html,
        Err(e) => {
            eprintln!("Template rendering error: {}", e);
            return;
        }
    };

    if let Some(path) = out_path {
        let mut file = File::create(path).expect("Could not create HTML output");
        file.write_all(rendered.as_bytes()).expect("Write failed");
    } else {
        println!("{}", rendered);
    }
}
