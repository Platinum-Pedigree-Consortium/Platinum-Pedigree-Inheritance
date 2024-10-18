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
    /// mothers name (matched to json)
    #[arg(short, long)]
    mother: String,
    /// father name (matched to json)
    #[arg(short, long)]
    father: String,
}

#[derive(Deserialize, Debug)]

struct Assembly {
    hap1: String,
    hap2: String,
}

struct Region {
    seqid: String,
    start: i64,
    end: i64,
}
impl Region {
    fn parse_region(region_str: &String) -> Result<Region, Box<dyn Error>> {
        // Parse the region
        let region_parts: Vec<&str> = region_str.split(':').collect();
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

        Ok((Region {
            seqid: chrom.to_string(),
            start: t_start,
            end: t_end,
        }))
    }
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
    region: &Region,
    flag: u16,
    flip_rc: bool,
) -> Result<Vec<FastaEntry>, Box<dyn Error>> {
    let mut bam = IndexedReader::from_path(in_bam)?;

    // Get reference sequence ID for chromosome name
    let tid = bam
        .header()
        .tid(region.seqid.as_bytes())
        .ok_or_else(|| format!("Chromosome '{}' not found in BAM file", region.seqid))?;

    // Fetch the reads that overlap with the region
    bam.fetch((tid, region.start, region.end))?;

    let mut fasta_entries: Vec<FastaEntry> = Vec::new();

    for r in bam.records() {
        let record = r.unwrap();

        // Map target region to query region (read coordinates)
        let cigar = record.cigar();
        let qpos = target_to_query(&cigar, record.pos(), region.start, region.end);

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

    Ok(fasta_entries)
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

fn test_concordance(
    haps: &HashMap<String, (String, String)>,
    iht: &InheritanceBlock,
    mother: &String,
    father: &String,
) -> bool {
    for i in 1..5 {
        println!("{} ", i);
        let mut lookup: HashMap<char, String> = HashMap::new();
        match i {
            1 => {
                lookup.insert('A', haps.get(father).unwrap().0.clone());
                lookup.insert('B', haps.get(father).unwrap().1.clone());
                lookup.insert('C', haps.get(mother).unwrap().0.clone());
                lookup.insert('D', haps.get(mother).unwrap().1.clone());
            }
            2 => {
                lookup.insert('A', haps.get(father).unwrap().1.clone());
                lookup.insert('B', haps.get(father).unwrap().0.clone());
                lookup.insert('C', haps.get(mother).unwrap().0.clone());
                lookup.insert('D', haps.get(mother).unwrap().1.clone());
            }
            3 => {
                lookup.insert('A', haps.get(father).unwrap().1.clone());
                lookup.insert('B', haps.get(father).unwrap().0.clone());
                lookup.insert('C', haps.get(mother).unwrap().1.clone());
                lookup.insert('D', haps.get(mother).unwrap().0.clone());
            }
            4 => {
                lookup.insert('A', haps.get(father).unwrap().0.clone());
                lookup.insert('B', haps.get(father).unwrap().1.clone());
                lookup.insert('C', haps.get(mother).unwrap().1.clone());
                lookup.insert('D', haps.get(mother).unwrap().0.clone());
            }
            _ => {}
        }

        let mut passing = 0;
        let mut total_count = 0;

        for (si, sample) in iht.samples.iter().enumerate() {
            let a = iht.parental_hap.get(si).unwrap().chars().nth(0).unwrap();
            let b = iht.parental_hap.get(si).unwrap().chars().nth(1).unwrap();

            let allele1 = lookup.get(&a).unwrap();
            let allele2 = lookup.get(&b).unwrap();

            let mut expected = (allele1, allele2);
            if allele1 > allele2 {
                expected.0 = allele2;
                expected.1 = allele1;
            }

            let sample_name = iht.samples.get(si).unwrap();

            if !haps.contains_key(sample_name) {
                continue;
            }
            total_count += 1;

            let seen = haps.get(sample_name).unwrap();
            if seen.0 == *expected.0 && seen.1 == *expected.1 {
                passing += 1;
            }

            println!(
                "{:?} {:?} {} {} e:{:?} s:{:?} p:{}",
                sample_name,
                iht.parental_hap.get(si).unwrap(),
                a,
                b,
                expected,
                seen,
                passing
            );
        }

        println!("n:{} passing:{}", total_count, passing);
    }

    return false;
}

fn main() {
    let flag = 260;
    let flip_rc = false;
    let args = Args::parse();
    let region = Region::parse_region(&args.region).unwrap();

    // Load up the inheritance vectors
    let mut inheritance = iht::parse_inht(args.inheritance);
    let mut current_block_idx: usize = 0;

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

        if h1_str < h2_str {
            loaded_haps.insert(sample.0.clone(), (h1_str, h2_str));
        } else {
            loaded_haps.insert(sample.0.clone(), (h2_str, h1_str));
        }
    }

    let mut block = iht::get_iht_block(
        &mut inheritance,
        &region.seqid,
        region.start.try_into().unwrap(),
        &mut current_block_idx,
    );

    // println!("\n{}\n", block.unwrap());

    test_concordance(&loaded_haps, block.unwrap(), &args.mother, &args.father);

    for person in loaded_haps {
        println!("{} {:?}", person.0, person.1);
    }
}
