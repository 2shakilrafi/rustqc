use serde::Serialize;
use std::collections::HashMap;

#[derive(Serialize)]
pub struct Report {
    pub filename: String,
    pub total_reads: usize,
    pub avg_length: f64,
    pub min_length: usize,
    pub max_length: usize,
    pub mode_length: usize,
    pub gc_content: f64,
    pub qualities_per_position: Vec<f64>,

    pub read_length_histogram: HashMap<usize, usize>,
    pub gc_percent_histogram: HashMap<u8, usize>,

    pub percent_a: Vec<f64>,
    pub percent_t: Vec<f64>,
    pub percent_g: Vec<f64>,
    pub percent_c: Vec<f64>,

    // ✅ NEW: avg quality per read histogram
    pub per_sequence_quality_histogram: HashMap<u8, usize>,
    pub n_content: Vec<f64>, // ← NEW

}
