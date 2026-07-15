use serde_json::Value;
use std::fs;
use std::path::PathBuf;
use std::process::Command;

#[test]
fn analyzes_fastq_and_writes_structured_report() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let input = manifest.join("tests/fixtures/tiny.fastq");
    let output = std::env::temp_dir().join(format!("rustqc-cli-{}.json", std::process::id()));

    let status = Command::new(env!("CARGO_BIN_EXE_rustqc"))
        .args([
            "--quiet",
            "--format",
            "json",
            "--phred-offset",
            "33",
            "--output",
        ])
        .arg(&output)
        .arg(&input)
        .status()
        .expect("rustqc should launch");
    assert!(status.success());

    let report: Value = serde_json::from_slice(&fs::read(&output).expect("report should exist"))
        .expect("report should be valid JSON");
    assert_eq!(report["basic_statistics"]["data"]["total_sequences"], 4);
    assert_eq!(
        report["basic_statistics"]["data"]["sequences_flagged_as_poor_quality"],
        1
    );
    assert_eq!(
        report["per_sequence_quality_scores"]["data"]["most_frequent_score"],
        40
    );
    assert_eq!(report["per_base_sequence_quality"]["data"][0]["mean"], 40.0);

    let _ = fs::remove_file(output);
}
