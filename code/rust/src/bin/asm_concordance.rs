use clap::Parser;
use rust_htslib::bam::record::{Cigar, CigarStringView};
use rust_htslib::bam::{IndexedReader, Read, Record};
use serde::Deserialize;
use serde_json::Value;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::error::Error;
use std::fs;

use concordance::iht;
use concordance::iht::InheritanceBlock;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]

struct Args {
    /// A formatted region string, e.g chr1:0-1000
    #[arg(short, long)]
    region: String,
    /// A json file containing a path to the aligned bam files
    #[arg(short, long)]
    samples: String,
    /// CSV containing the inheritance vectors (optionally gzipped)
    #[arg(short, long)]
    inheritance: String,
}
#[derive(Deserialize, Debug)]

struct Assembly {
    hap1: String,
    hap2: String,
}

fn target_to_query(cigar: &CigarStringView, t_pos: i64, t_start: i64, t_end: i64) -> (i64, i64) {
    let mut t_off = t_pos;
    let mut q_off: i64 = 0;
    let mut q_start: i64 = -1;
    let mut q_end: i64 = -1;
    let mut start_unset = true;
    let mut end_unset = true;

    for op in cigar.iter() {
        let len = op.len() as i64;
        match *op {
            // Match or Mismatch (both query and reference advance)
            Cigar::Match(_) | Cigar::Equal(_) | Cigar::Diff(_) => {
                if t_off + len > t_start && start_unset {
                    q_start = q_off + (t_start - t_off).max(0);
                    start_unset = false;
                }
                if t_off + len >= t_end && end_unset {
                    q_end = q_off + (t_end - t_off).min(len);
                    end_unset = false;
                }
                q_off += len;
                t_off += len;
            }

            // Insertion (query advances, reference doesn't)
            Cigar::Ins(_) => {
                q_off += len;
            }

            // Deletion (reference advances, query doesn't)
            Cigar::Del(_) => {
                if t_off + len > t_start && start_unset {
                    q_start = q_off;
                    start_unset = false;
                }
                if t_off + len >= t_end && end_unset {
                    q_end = q_off;
                    end_unset = false;
                }
                t_off += len;
            }

            // Soft Clipping (query advances, reference doesn't)
            Cigar::SoftClip(_) => {
                q_off += len;
            }

            // Hard Clipping (query advances, reference doesn't, but it doesn't consume reference)
            Cigar::HardClip(_) => {
                // Hard Clipping does not affect query position for coordinate calculation
            }

            // Skip (advances reference, but not query)
            Cigar::RefSkip(_) => {
                t_off += len;
            }

            _ => {
                // Handle any other CIGAR operations if necessary
            }
        }

        // Exit early if both start and end positions are found
        if !start_unset && !end_unset {
            break;
        }
    }

    // If q_end is still unset, set it to q_off (end of query sequence)
    if q_end == -1 {
        q_end = q_off;
    }

    (q_start, q_end)
}

// Struct for holding extracted fasta entries
#[derive(Debug)]
struct FastaEntry {
    strand: bool,
    length: usize,
    name: String,
    seq: String,
}

// Compare FastaEntries by strand
fn my_cmp(x: &FastaEntry, y: &FastaEntry) -> Ordering {
    x.strand.cmp(&y.strand)
}

// Function to extract sequences from BAM, trim to target region, and optionally reverse complement
fn bam_to_fasta(
    in_bam: &str,
    region: &str,
    flag: u16,
    flip_rc: bool,
) -> Result<(Vec<FastaEntry>), Box<dyn Error>> {
    let mut bam = IndexedReader::from_path(in_bam)?;

    // Parse the region
    let region_parts: Vec<&str> = region.split(':').collect();
    if region_parts.len() != 2 {
        return Err("Region format is incorrect, expected 'chr:start-end'".into());
    }

    let chrom = region_parts[0];
    let pos_parts: Vec<&str> = region_parts[1].split('-').collect();
    if pos_parts.len() != 2 {
        return Err("Region format is incorrect, expected 'start-end'".into());
    }

    let t_start: i64 = pos_parts[0].parse()?;
    let t_end: i64 = pos_parts[1].parse()?;

    // Get reference sequence ID for chromosome name
    let tid = bam
        .header()
        .tid(chrom.as_bytes())
        .ok_or_else(|| format!("Chromosome '{}' not found in BAM file", chrom))?;

    // Fetch the reads that overlap with the region
    bam.fetch((tid, t_start, t_end))?;

    let mut fasta_entries: Vec<FastaEntry> = Vec::new();

    for r in bam.records() {
        let record = r.unwrap();

        // Map target region to query region (read coordinates)
        let cigar = record.cigar();
        let qpos = target_to_query(&cigar, record.pos(), t_start, t_end);

        // Extract the subsequence within the query region
        let seq = record.seq().as_bytes();
        let q_start = qpos.0 as usize;
        let q_end = qpos.1 as usize;

        if q_start >= seq.len() || q_end >= seq.len() || q_start > q_end {
            continue;
        }

        let mut dna_seq = String::from_utf8(seq[q_start..=q_end].to_vec())?;

        // Optionally reverse complement the sequence
        if flip_rc && record.is_reverse() {
            dna_seq = reverse_complement(&dna_seq);
        }

        // Add the fasta entry
        fasta_entries.push(FastaEntry {
            strand: record.is_reverse(),
            length: dna_seq.len(),
            name: String::from_utf8(record.qname().to_vec())?,
            seq: dna_seq,
        });
    }

    // Sort fasta entries by strand
    fasta_entries.sort_by(my_cmp);

    Ok((fasta_entries))
}

// Reverse complement function for DNA sequences
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => 'N',
        })
        .collect()
}

fn main() {
    let flag = 260;
    let flip_rc = false;
    let args = Args::parse();
    let region = args.region;

    // Load up the inheritance vectors
    let mut inheritance = iht::parse_inht(args.inheritance);

    // Read the JSON file into a string
    let json_data = fs::read_to_string(args.samples).unwrap();

    // Parse the JSON string
    let sample_info: HashMap<String, Assembly> = serde_json::from_str(&json_data).unwrap();

    let mut loaded_haps: HashMap<String, (String, String)> = HashMap::new();

    for sample in sample_info {
        let h1 = bam_to_fasta(&sample.1.hap1, &region, flag, flip_rc).unwrap();
        let h2 = bam_to_fasta(&sample.1.hap2, &region, flag, flip_rc).unwrap();

        let mut h1_str = "".to_string();
        if h1.len() == 1 {
            h1_str = h1.get(0).unwrap().seq.clone();
        }

        let mut h2_str = "".to_string();
        if h2.len() == 1 {
            h2_str = h2.get(0).unwrap().seq.clone();
        }

        loaded_haps.insert(sample.0.clone(), (h1_str, h2_str));
    }

    for person in loaded_haps {
        println!("{} {:?}", person.0, person.1);
    }
}
