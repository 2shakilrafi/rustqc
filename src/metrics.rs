use crate::report::Status;
use std::fs;
use std::path::Path;

pub const OVERREP_WARN_PERCENT: f64 = 0.1;
pub const OVERREP_FAIL_PERCENT: f64 = 1.0;
pub const MAX_TRACKED_UNIQUE_SEQUENCES: usize = 100_000;
pub const MAX_TILES: usize = 2_500;

#[derive(Clone, Debug)]
pub struct PositionGroup {
    pub start: usize,
    pub end: usize,
}

impl PositionGroup {
    pub fn label(&self) -> String {
        if self.start == self.end {
            self.start.to_string()
        } else {
            format!("{}-{}", self.start, self.end)
        }
    }

    pub fn midpoint(&self) -> f64 {
        (self.start + self.end) as f64 / 2.0
    }
}

pub fn make_base_groups(
    max_length: usize,
    no_group: bool,
    min_length: usize,
) -> Vec<PositionGroup> {
    let max_length = max_length.max(min_length);
    if max_length == 0 {
        return Vec::new();
    }
    if no_group || max_length <= 75 {
        return (1..=max_length)
            .map(|position| PositionGroup {
                start: position,
                end: position,
            })
            .collect();
    }

    let interval = linear_interval(max_length);
    let mut groups = Vec::with_capacity(75);
    let mut start = 1;
    while start <= max_length {
        let mut end = start + interval - 1;
        if start < 10 {
            end = start;
        } else if start == 10 && interval > 10 {
            end = interval - 1;
        }
        end = end.min(max_length);
        groups.push(PositionGroup { start, end });
        if start < 10 {
            start += 1;
        } else if start == 10 && interval > 10 {
            start = interval;
        } else {
            start += interval;
        }
    }
    groups
}

fn linear_interval(length: usize) -> usize {
    let mut multiplier = 1;
    loop {
        for base in [2, 5, 10] {
            let interval = base * multiplier;
            let tail = length.saturating_sub(9);
            let group_count = 9 + tail.div_ceil(interval);
            if group_count < 75 {
                return interval;
            }
        }
        multiplier *= 10;
    }
}

pub fn high_is_bad(value: f64, warn: f64, fail: f64) -> Status {
    if value > fail {
        Status::Fail
    } else if value > warn {
        Status::Warn
    } else {
        Status::Pass
    }
}

pub fn low_is_bad(value: f64, warn: f64, fail: f64) -> Status {
    if value < fail {
        Status::Fail
    } else if value < warn {
        Status::Warn
    } else {
        Status::Pass
    }
}

#[derive(Clone, Debug)]
pub struct NamedSequence {
    pub name: String,
    pub sequence: Vec<u8>,
}

pub fn load_named_sequences(path: &Path) -> Result<Vec<NamedSequence>, String> {
    let text = fs::read_to_string(path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    parse_named_sequences(&text)
}

pub fn parse_named_sequences(text: &str) -> Result<Vec<NamedSequence>, String> {
    let mut sequences = Vec::new();
    for (line_number, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some((name, sequence)) = line.split_once('\t') else {
            return Err(format!(
                "line {} must contain name<TAB>sequence",
                line_number + 1
            ));
        };
        let sequence = sequence.trim().as_bytes().to_ascii_uppercase();
        if sequence.is_empty()
            || !sequence
                .iter()
                .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
        {
            return Err(format!(
                "line {} contains an invalid DNA sequence",
                line_number + 1
            ));
        }
        sequences.push(NamedSequence {
            name: name.trim().to_string(),
            sequence,
        });
    }
    if sequences.is_empty() {
        return Err("sequence list is empty".into());
    }
    Ok(sequences)
}

pub fn default_adapters() -> Vec<NamedSequence> {
    parse_named_sequences(
        "Illumina Universal Adapter\tAGATCGGAAGAG\n\
         Illumina Small RNA 3' Adapter\tTGGAATTCTCGG\n\
         Illumina Small RNA 5' Adapter\tGATCGTCGGACT\n\
         Nextera Transposase Sequence\tCTGTCTCTTATA\n\
         PolyA\tAAAAAAAAAAAA\n\
         PolyG\tGGGGGGGGGGGG\n",
    )
    .expect("built-in adapter definitions are valid")
}

pub fn default_contaminants() -> Vec<NamedSequence> {
    parse_named_sequences(
        "Illumina Single End Adapter 1\tGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG\n\
         Illumina Single End Adapter 2\tCAAGCAGAAGACGGCATACGAGCTCTTCCGATCT\n\
         Illumina Paired End Adapter 1\tACACTCTTTCCCTACACGACGCTCTTCCGATCT\n\
         Illumina Paired End Adapter 2\tGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG\n\
         Illumina Small RNA Adapter 1\tGTTCAGAGTTCTACAGTCCGACGATC\n\
         Illumina Small RNA Adapter 2\tTGGAATTCTCGGGTGCCAAGG\n\
         Illumina Multiplexing Adapter 1\tGATCGGAAGAGCACACGTCT\n\
         Illumina Multiplexing Adapter 2\tACACTCTTTCCCTACACGACGCTCTTCCGATCT\n\
         Illumina Multiplexing PCR Primer 2.01\tGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT\n\
         TruSeq Universal Adapter\tAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n\
         TruSeq Adapter, Index 1\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 2\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 3\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 4\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 5\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 6\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 7\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 8\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 9\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 10\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 11\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG\n\
         TruSeq Adapter, Index 12\tGATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG\n",
    )
    .expect("built-in contaminant definitions are valid")
}

#[derive(Clone, Copy, Debug)]
struct Match {
    length: usize,
    identity: usize,
}

pub fn possible_source(query: &[u8], contaminants: &[NamedSequence]) -> String {
    let query = query.to_ascii_uppercase();
    let mut best: Option<(&str, Match)> = None;
    for contaminant in contaminants {
        let reverse = reverse_complement(&contaminant.sequence);
        for target in [&contaminant.sequence[..], &reverse[..]] {
            let Some(candidate) = best_match(&query, target) else {
                continue;
            };
            let replace = best.is_none_or(|(_, current)| {
                candidate.length > current.length
                    || (candidate.length == current.length && candidate.identity > current.identity)
            });
            if replace {
                best = Some((&contaminant.name, candidate));
            }
        }
    }
    best.map_or_else(
        || "No Hit".to_string(),
        |(name, hit)| format!("{name} ({}% over {}bp)", hit.identity, hit.length),
    )
}

fn best_match(query: &[u8], target: &[u8]) -> Option<Match> {
    if (8..20).contains(&query.len()) && target.windows(query.len()).any(|window| window == query) {
        return Some(Match {
            length: query.len(),
            identity: 100,
        });
    }
    if query.len() < 20 || target.len() < 20 {
        return None;
    }

    let mut best = None;
    let min_offset = -(target.len() as isize - 20);
    let max_offset = query.len() as isize - 20;
    for offset in min_offset..=max_offset {
        let target_start = (-offset).max(0) as usize;
        let query_start = offset.max(0) as usize;
        let overlap = (target.len() - target_start).min(query.len() - query_start);
        if overlap < 20 {
            continue;
        }

        let mut left = 0;
        let mut mismatches = 0;
        for right in 0..overlap {
            if target[target_start + right] != query[query_start + right] {
                mismatches += 1;
            }
            while mismatches > 1 {
                if target[target_start + left] != query[query_start + left] {
                    mismatches -= 1;
                }
                left += 1;
            }
            let length = right + 1 - left;
            if length >= 20 {
                let identity = ((length - mismatches) * 100) / length;
                let candidate = Match { length, identity };
                if best.is_none_or(|current: Match| {
                    candidate.length > current.length
                        || (candidate.length == current.length
                            && candidate.identity > current.identity)
                }) {
                    best = Some(candidate);
                }
            }
        }
    }
    best
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups_150_base_reads_like_fastqc() {
        let groups = make_base_groups(150, false, 0);
        let labels: Vec<_> = groups.iter().map(PositionGroup::label).collect();
        assert_eq!(
            &labels[..10],
            ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10-14"]
        );
        assert_eq!(labels.last().unwrap(), "150");
        assert_eq!(labels.len(), 38);
    }

    #[test]
    fn identifies_truseq_index() {
        let source = possible_source(
            b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGC",
            &default_contaminants(),
        );
        assert!(source.starts_with("TruSeq Adapter, Index 8 (100% over 50bp)"));
    }
}
