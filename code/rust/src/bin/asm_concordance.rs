use clap::Parser;
use log::debug;
use log::info;
use log::warn;
use log::LevelFilter;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Cigar, CigarStringView};
use rust_htslib::bam::{IndexedReader, Read, Record};
use serde::Deserialize;
use serde_json::Value;
use std::cmp::Ordering;
use std::collections::{hash_map, HashMap};
use std::convert::TryInto;
use std::error;
use std::error::Error;
use std::fmt;
use std::fs;

use concordance::bed;
use concordance::bed::BedRecord;
use concordance::iht;
use concordance::iht::InheritanceBlock;

type Result<T> = std::result::Result<T, Box<dyn error::Error>>;

use rust_lapper::{Interval, Lapper};

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
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
}

#[derive(Deserialize, Debug)]

struct Assembly {
    hap1: String,
    hap2: String,
    sex: String,
}

#[derive(Clone)]
struct record_bucket {
    rec: Record,
    seq: Vec<u8>,
}

impl PartialEq for record_bucket {
    fn eq(&self, other: &Self) -> bool {
        self.rec == other.rec
    }
}

impl Eq for record_bucket {}

#[derive(Deserialize, Debug)]

struct Region {
    seqid: String,
    start: i64,
    end: i64,
}

impl From<BedRecord> for Region {
    fn from(bed_record: BedRecord) -> Self {
        Region {
            seqid: bed_record.chrom,
            start: bed_record.start as i64, // Convert u64 to i64
            end: bed_record.end as i64,     // Convert u64 to i64
        }
    }
}

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}\t{}\t{}", self.seqid, self.start, self.end)
    }
}

impl Region {
    fn parse_region(region_str: &String) -> Result<Region> {
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
#[derive(Debug, Clone)]
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
fn data_to_fasta(
    data: &mut Lapper<usize, record_bucket>,
    region: &Region,
    flag: u16,
    flip_rc: bool,
) -> Result<Vec<FastaEntry>> {
    let mut fasta_entries: Vec<FastaEntry> = Vec::new();

    for r in data.find(region.start as usize, region.end as usize) {
        let record = &r.val.rec;

        // Map target region to query region (read coordinates)
        let cigar = record.cigar();
        let qpos = target_to_query(&cigar, record.pos(), region.start, region.end);

        // Extract the subsequence within the query region

        let q_start = qpos.0 as usize;
        let q_end = qpos.1 as usize;

        if q_start >= r.val.seq.len() || q_end >= r.val.seq.len() || q_start > q_end {
            continue;
        }

        let mut dna_seq = String::from_utf8(r.val.seq[q_start..=q_end].to_vec())?;

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

    return Ok(fasta_entries);
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

        // A,B C,D
        // AA,TC AA,TG

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
            /*
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
            */
        }

        if passing == total_count && total_count != 0 {
            return true;
        }

        // println!("n:{} passing:{}", total_count, passing);
    }

    return false;
}

fn open_bam_files(fns: &HashMap<String, Assembly>) -> Result<HashMap<String, [IndexedReader; 2]>> {
    let mut open_bh: HashMap<String, [IndexedReader; 2]> = HashMap::new();

    for s in fns {
        open_bh.insert(
            s.0.clone(),
            [
                IndexedReader::from_path(s.1.hap1.clone()).unwrap(),
                IndexedReader::from_path(s.1.hap2.clone()).unwrap(),
            ],
        );
    }
    return Ok(open_bh);
}

fn fetch_haplotypes(
    data: &mut HashMap<String, [Lapper<usize, record_bucket>; 2]>,
    region: &Region,
    flag: u16,
    flip_rc: bool,
    problem_counts: &mut (i32, i32),
    sample_info: &HashMap<String, Assembly>,
) -> Result<HashMap<String, (String, String)>> {
    let mut loaded_haps: HashMap<String, (String, String)> = HashMap::new();

    for sample in data {
        let sex = sample_info.get(sample.0).unwrap().sex.clone();

        let mut h1 = data_to_fasta(sample.1.get_mut(0).unwrap(), &region, flag, flip_rc).unwrap();
        let mut h2 = data_to_fasta(sample.1.get_mut(1).unwrap(), &region, flag, flip_rc).unwrap();

        if sex == "male"
            && (region.seqid == "X" || region.seqid == "chrX" || region.seqid == "ChrX")
        {
            if h1.len() == 0 && h2.len() == 1 {
                h1 = h2.to_vec();
            }
            if h2.len() == 0 && h1.len() == 1 {
                h2 = h1.to_vec();
            }
        }

        if h1.len() > 1 {
            problem_counts.1 += 1;
            continue;
        }
        if h2.len() > 1 {
            problem_counts.1 += 1;
            continue;
        }

        if h1.len() < 1 {
            problem_counts.0 += 1;
            continue;
        }
        if h2.len() < 1 {
            problem_counts.0 += 1;
            continue;
        }

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
    Ok(loaded_haps)
}

fn load_chromosome(
    samples: &mut HashMap<String, [IndexedReader; 2]>,
    chr: &String,
    data: &mut HashMap<String, [Lapper<usize, record_bucket>; 2]>,
) {
    data.clear();

    type Iv = Interval<usize, record_bucket>;
    let mut storage: HashMap<String, [Vec<Iv>; 2]> = HashMap::new();

    // Initialize the storage HashMap with empty Vecs for each sample
    for sample in samples.iter() {
        storage.insert(sample.0.clone(), [Vec::new(), Vec::new()]);
    }

    // Iterate over the samples and fetch records for the given chromosome
    for sample in samples.iter_mut() {
        for bi in 0..2 {
            // Get a mutable reference to the BAM reader using `get_mut`
            let bam: &mut IndexedReader = sample.1.get_mut(bi).unwrap();

            // Find the target ID (tid) for the chromosome in the BAM header
            let tid = bam
                .header()
                .tid(chr.as_bytes())
                .ok_or_else(|| format!("Chromosome '{}' not found in BAM file", chr))
                .unwrap();

            // Fetch the reads for the specific chromosome
            bam.fetch(tid).unwrap();

            info!(
                "loading mapped contigs for sample: {} haplotype: {} chromosome: {}",
                sample.0, bi, chr
            );
            // Iterate over the records in the BAM file and push them to storage
            for r in bam.records() {
                let rec = r.unwrap();
                let seq = rec.seq().as_bytes();

                let val = record_bucket { rec, seq };
                let s = storage.get_mut(&sample.0.clone()).unwrap();

                let start = val.rec.reference_start() as usize;
                let stop = val.rec.reference_end() as usize;
                s[bi].push(Iv { start, stop, val });
            }
        }
    }
    for sample in storage {
        data.insert(
            sample.0,
            [
                Lapper::new((*sample.1.get(0).unwrap().clone()).to_vec()),
                Lapper::new((*sample.1.get(1).unwrap().clone()).to_vec()),
            ],
        );
    }
}
fn main() {
    let flag = 260;
    let flip_rc = false;
    let args = Args::parse();

    let filter_level: LevelFilter = match args.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    // Load up the inheritance vectors
    let mut inheritance = iht::parse_inht(args.inheritance);

    // Read the JSON file into a string
    let json_data = fs::read_to_string(args.samples).unwrap();

    // Parse the JSON string
    let sample_info: HashMap<String, Assembly> = serde_json::from_str(&json_data).unwrap();

    // open bam files
    let mut bams: HashMap<String, [IndexedReader; 2]> = open_bam_files(&sample_info).unwrap();

    let bed_records = concordance::bed::read_bed_file(args.region).unwrap();

    let mut last_chrom = "".to_string();

    let mut read_data: HashMap<String, [Lapper<usize, record_bucket>; 2]> = HashMap::new();

    for region in bed_records {
        let lr: Region = region.into();
        if last_chrom != lr.seqid {
            last_chrom = lr.seqid.clone();
            load_chromosome(&mut bams, &last_chrom, &mut read_data);
        }

        let mut problems = (0, 0);
        let loaded_haps = fetch_haplotypes(
            &mut read_data,
            &lr,
            260 as u16,
            false,
            &mut problems,
            &sample_info,
        );

        let mut current_block_idx: usize = 0;
        let mut block = iht::get_iht_block(
            &mut inheritance,
            &lr.seqid,
            lr.start.try_into().unwrap(),
            &mut current_block_idx,
        );

        let mut failed_block = false;
        match block {
            None => {
                warn!(
                    "Warning: skipping variant as it is not in a block {} {}",
                    lr.seqid, lr.start
                );
                failed_block = true;
            }
            Some(_) => {}
        }

        if problems.0 > 0 || problems.1 > 0 || failed_block {
            println!(
                "{}\t{}\t{}\t{}\t{}",
                lr, false, problems.0, problems.1, failed_block
            );
            continue;
        }

        let concordant = test_concordance(
            &loaded_haps.unwrap(),
            block.unwrap(),
            &args.mother,
            &args.father,
        );

        println!("{}\t{}\t0\t0\t{}", lr, concordant, failed_block);
    }
}
