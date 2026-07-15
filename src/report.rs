use serde::Serialize;
use std::fmt;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Status {
    Pass,
    Warn,
    Fail,
    Skip,
}

impl fmt::Display for Status {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Pass => "PASS",
            Self::Warn => "WARN",
            Self::Fail => "FAIL",
            Self::Skip => "SKIP",
        })
    }
}

#[derive(Debug, Serialize)]
pub struct Module<T> {
    pub status: Status,
    pub data: T,
}

#[derive(Debug, Serialize)]
pub struct Report {
    pub schema_version: u8,
    pub program: &'static str,
    pub version: &'static str,
    pub generated_at_unix: u64,
    pub basic_statistics: Module<BasicStatistics>,
    pub per_base_sequence_quality: Module<Vec<BaseQualityPoint>>,
    pub per_tile_sequence_quality: Module<TileQualityData>,
    pub per_sequence_quality_scores: Module<SequenceQualityData>,
    pub per_base_sequence_content: Module<Vec<BaseContentPoint>>,
    pub per_sequence_gc_content: Module<GcContentData>,
    pub per_base_n_content: Module<Vec<PositionValue>>,
    pub sequence_length_distribution: Module<Vec<CountPoint>>,
    pub sequence_duplication_levels: Module<DuplicationData>,
    pub overrepresented_sequences: Module<Vec<OverrepresentedSequence>>,
    pub adapter_content: Module<AdapterContentData>,
}

impl Report {
    pub fn summaries(&self) -> [(&'static str, Status); 11] {
        [
            ("Basic Statistics", self.basic_statistics.status),
            (
                "Per base sequence quality",
                self.per_base_sequence_quality.status,
            ),
            (
                "Per tile sequence quality",
                self.per_tile_sequence_quality.status,
            ),
            (
                "Per sequence quality scores",
                self.per_sequence_quality_scores.status,
            ),
            (
                "Per base sequence content",
                self.per_base_sequence_content.status,
            ),
            (
                "Per sequence GC content",
                self.per_sequence_gc_content.status,
            ),
            ("Per base N content", self.per_base_n_content.status),
            (
                "Sequence Length Distribution",
                self.sequence_length_distribution.status,
            ),
            (
                "Sequence Duplication Levels",
                self.sequence_duplication_levels.status,
            ),
            (
                "Overrepresented sequences",
                self.overrepresented_sequences.status,
            ),
            ("Adapter Content", self.adapter_content.status),
        ]
    }
}

#[derive(Debug, Serialize)]
pub struct BasicStatistics {
    pub filename: String,
    pub file_type: &'static str,
    pub encoding: String,
    pub total_sequences: u64,
    pub total_bases: u64,
    pub sequences_flagged_as_poor_quality: u64,
    pub sequence_length: String,
    pub min_length: usize,
    pub max_length: usize,
    pub mean_length: f64,
    pub gc_percent: f64,
}

#[derive(Debug, Serialize)]
pub struct BaseQualityPoint {
    pub base: String,
    pub position: f64,
    pub mean: f64,
    pub median: Option<f64>,
    pub lower_quartile: Option<f64>,
    pub upper_quartile: Option<f64>,
    pub percentile_10: Option<f64>,
    pub percentile_90: Option<f64>,
}

#[derive(Debug, Default, Serialize)]
pub struct TileQualityData {
    pub tiles: Vec<u32>,
    pub positions: Vec<String>,
    pub position_midpoints: Vec<f64>,
    pub deviations: Vec<Vec<f64>>,
    pub max_deviation: f64,
}

#[derive(Debug, Serialize)]
pub struct SequenceQualityData {
    pub distribution: Vec<CountPoint>,
    pub most_frequent_score: i32,
}

#[derive(Debug, Serialize)]
pub struct BaseContentPoint {
    pub base: String,
    pub position: f64,
    pub g: f64,
    pub a: f64,
    pub t: f64,
    pub c: f64,
}

#[derive(Debug, Serialize)]
pub struct GcContentData {
    pub distribution: Vec<GcPoint>,
    pub deviation_percent: f64,
}

#[derive(Debug, Serialize)]
pub struct GcPoint {
    pub gc_percent: u8,
    pub count: f64,
    pub theoretical: f64,
}

#[derive(Debug, Serialize)]
pub struct PositionValue {
    pub base: String,
    pub position: f64,
    pub value: f64,
}

#[derive(Debug, Serialize)]
pub struct CountPoint {
    pub value: u64,
    pub label: String,
    pub count: f64,
}

#[derive(Debug, Serialize)]
pub struct DuplicationData {
    pub total_deduplicated_percent: f64,
    pub levels: Vec<DuplicationPoint>,
}

#[derive(Debug, Serialize)]
pub struct DuplicationPoint {
    pub level: String,
    pub percentage_of_total: f64,
}

#[derive(Debug, Serialize)]
pub struct OverrepresentedSequence {
    pub sequence: String,
    pub count: u64,
    pub percentage: f64,
    pub possible_source: String,
}

#[derive(Debug, Serialize)]
pub struct AdapterContentData {
    pub positions: Vec<String>,
    pub position_midpoints: Vec<f64>,
    pub series: Vec<AdapterSeries>,
    pub max_content: f64,
}

#[derive(Debug, Serialize)]
pub struct AdapterSeries {
    pub name: String,
    pub percentages: Vec<f64>,
}
