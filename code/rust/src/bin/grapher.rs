use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::read_to_string;
use std::{fmt, fs::File, io::BufWriter, io::Write, str};

use clap::Parser;

use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::{record::GenotypeAllele, IndexedReader, Read};
use rust_htslib::faidx::Reader as FaidxReader;

use log::info;
use serde::Serialize;

/// A tool for printing all possible haplotypes
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// fasta
    #[arg(short, long)]
    fasta: String,

    /// VCF file
    #[arg(short, long)]
    vcf: String,

    /// sample name (must be in vcf header)
    #[arg(short, long)]
    sample: String,

    /// region (e.g. chr1:1-500)
    #[arg(short, long)]
    region: String,

    /// prefix for output
    #[arg(short, long)]
    prefix: String,

    /// discover only
    #[arg(short, long, action)]
    discovery_only: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
enum VarType {
    Snv,
    Insertion,
    Deletion,
    Ref,
}

#[derive(Debug, PartialEq, Eq, Clone, Serialize)]
enum OvlType {
    NoOvl,
    OvlNext,
    OvlPrev,
}
#[derive(Clone)]
struct Region {
    rid: u32,
    name: String,
    start: u64,
    end: u64,
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.name, self.start, self.end)
    }
}

struct Data {
    region: Region,
    variants: [Vec<Var>; 2],
    sequence: String,
}

#[derive(Debug, Clone, Serialize)]
struct Var {
    seqid: String,
    vt: VarType,
    start: i64,
    end: i64,
    ref_allele: String,
    alt_allele: String,
    ovl: OvlType,
    next_vars: Vec<usize>,
    idx: usize,
    vidx: usize,
}

#[derive(Debug, Clone, Serialize)]
struct VarTainer {
    label: String,
    count: usize,
    vars: Vec<Var>,
}

impl fmt::Display for Var {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{:?}\t{}\t{}",
            self.seqid,
            self.start + 1,
            self.ref_allele,
            self.alt_allele,
            self.ovl,
            self.idx,
            self.vidx,
        )
    }
}

impl Var {
    fn varlen(&self) -> usize {
        return self.ref_allele.len() - self.alt_allele.len();
    }
}

fn get_end(pos: i64, vtype: VarType, alleles: Vec<&[u8]>, alt_idx: usize) -> i64 {
    assert!(alleles.len() > alt_idx);
    if vtype == VarType::Snv || vtype == VarType::Insertion {
        return pos;
    }
    return pos + (alleles[0].len() - alleles[alt_idx].len()) as i64;
}

fn get_var_type(alleles: Vec<&[u8]>, alt_idx: usize) -> VarType {
    assert!(alleles.len() > alt_idx);
    if alleles[0].len() == alleles[alt_idx].len() && alleles[alt_idx].len() == 1 {
        return VarType::Snv;
    } else if alleles[0].len() > alleles[alt_idx].len() {
        return VarType::Deletion;
    } else {
        return VarType::Insertion;
    }
}

fn parse_region(rstring: String, header: &HeaderView) -> Region {
    let region_parts: Vec<String> = rstring.split([':', '-']).map(str::to_string).collect();

    Region {
        rid: header.name2rid(region_parts[0].as_bytes()).unwrap(),
        name: region_parts[0].clone(),
        start: region_parts[1].parse::<u64>().unwrap(),
        end: region_parts[2].parse::<u64>().unwrap(),
    }
}

fn build_haplotypes(haplotype: &Vec<Var>, refseq: &String, region_start: i64) -> String {
    let mut hap: Vec<char> = refseq.clone().chars().collect();
    for v in haplotype {
        match v.vt {
            VarType::Snv => {
                let ref_base = v.ref_allele.as_bytes()[0] as char;
                let seq_base = hap[(v.start - region_start as i64) as usize].to_ascii_uppercase();
                assert!(ref_base == seq_base);
                hap[(v.start - region_start as i64) as usize] = v.alt_allele.as_bytes()[0] as char
            }
            VarType::Deletion => {
                let left = (v.start + 1) - region_start as i64;
                let mut right = left + v.varlen() as i64;

                right = std::cmp::min(right, (hap.len() as i64) - 1);

                for position in left..right {
                    hap[position as usize] = '-';
                }
            }
            // this will catch REF types
            _ => {}
        }
    }
    for i in haplotype.iter().rev() {
        if i.vt == VarType::Insertion {
            let mut insert_seq = i.alt_allele.chars().collect::<Vec<char>>();
            insert_seq.remove(0);
            let left = (i.start + 1) - region_start as i64;
            hap.splice((left as usize)..(left as usize), insert_seq);
        }
    }
    hap.iter().filter(|dna_base| **dna_base != '-').collect()
}

fn load_data(
    args: &Args,
    bcf: &mut IndexedReader,
    region: &Region,
    sample_lookup: &HashMap<String, usize>,
) -> Data {
    let mut variants = [Vec::new(), Vec::new()];

    let fasta = FaidxReader::from_path((*args).fasta.clone()).expect("Failed to open FASTA");

    let region_seq = fasta
        .fetch_seq_string(
            region.name.clone(),
            region.start as usize,
            region.end as usize,
        )
        .expect("FAILED to get fasta region")
        .to_ascii_lowercase();

    let region_bytes = region_seq.as_bytes();

    bcf.fetch(region.rid, region.start, Some(region.end))
        .unwrap();

    let mut dummy = Var {
        seqid: region.name.clone(),
        vt: VarType::Ref,
        start: region.start as i64,
        end: region.start as i64,
        ref_allele: "".to_string(),
        alt_allele: "".to_string(),
        ovl: OvlType::NoOvl,
        next_vars: Vec::new(),
        idx: 0 as usize,
        vidx: 0 as usize,
    };

    variants[0].push(dummy.clone());
    dummy.vidx = 1;
    variants[1].push(dummy.clone());

    for record in bcf.records() {
        let r = record.unwrap();

        let gt = r
            .genotypes()
            .unwrap()
            .get(*sample_lookup.get(&args.sample).unwrap());

        for vindex in 0..2 {
            // phased or unphased reference alleles are skipped, as they do not mutate reference fasta
            if gt[vindex] == GenotypeAllele::Unphased(0) || gt[vindex] == GenotypeAllele::Phased(0)
            {
                continue;
            }

            let ref_idx = 0 as usize;
            let alt_idx = gt[vindex].index().unwrap() as usize;

            let vtype = get_var_type(r.alleles(), alt_idx);
            let end = get_end(r.pos(), vtype.clone(), r.alleles(), alt_idx);

            let mut last_idx: usize = variants[vindex].len() - 1;
            let current_idx: usize = variants[vindex].len();

            let mut entry = Var {
                seqid: region.name.clone(),
                vt: vtype.clone(),
                start: r.pos(),
                end: end,
                ref_allele: (*str::from_utf8(r.alleles()[ref_idx]).unwrap()).to_string(),
                alt_allele: (*str::from_utf8(r.alleles()[alt_idx]).unwrap()).to_string(),
                ovl: OvlType::NoOvl,
                next_vars: Vec::new(),
                idx: current_idx,
                vidx: vindex,
            };

            if entry.vt == VarType::Deletion || entry.vt == VarType::Insertion {
                entry.start += 1;
            }

            // checking that everything lines up with the reference
            if entry.vt == VarType::Snv {
                let base = (region_bytes[(entry.start - (region.start as i64)) as usize] as char)
                    .to_ascii_uppercase();
                assert!(entry.ref_allele.as_bytes()[0] as char == base);
            }

            let mut last_var: &mut Var = &mut variants[vindex][last_idx];

            if last_var.end >= entry.start {
                (*last_var).ovl = OvlType::OvlNext;
                entry.ovl = OvlType::OvlPrev;
            }

            while true {
                if last_var.end < entry.start {
                    last_var.next_vars.push(current_idx);
                }

                if last_var.ovl == OvlType::NoOvl || last_idx == 0 {
                    break;
                }
                last_idx -= 1;
                last_var = &mut variants[vindex][last_idx];
            }

            if entry.vt == VarType::Deletion || entry.vt == VarType::Insertion {
                entry.start -= 1;
            }

            variants[vindex].push(entry);
        }
    }
    return Data {
        region: Region {
            rid: region.rid,
            name: region.name.clone(),
            start: region.start,
            end: region.end,
        },
        variants: variants,
        sequence: region_seq,
    };
}

fn find_all_paths(graph: &[Var], start_idx: usize) -> Vec<Vec<Var>> {
    let mut paths = Vec::new();
    let mut visited = HashSet::new();

    dfs(&graph, start_idx, &mut vec![], &mut visited, &mut paths);

    paths
}

fn dfs(
    graph: &[Var],
    current_idx: usize,
    current_path: &mut Vec<Var>,
    visited: &mut HashSet<usize>,
    all_paths: &mut Vec<Vec<Var>>,
) {
    visited.insert(current_idx);
    current_path.push(graph[current_idx].clone());

    if graph[current_idx].next_vars.is_empty() {
        // Reached a leaf node, save the current path
        all_paths.push(current_path.clone());
    } else {
        for &next_idx in &graph[current_idx].next_vars {
            if !visited.contains(&next_idx) {
                dfs(graph, next_idx, current_path, visited, all_paths);
            }
        }
    }

    // Backtrack
    visited.remove(&current_idx);
    current_path.pop();
}

fn parse_regions(args: &Args) -> Vec<Region> {
    let mut regions: Vec<Region> = Vec::new();
    let bcf = IndexedReader::from_path((*args).vcf.clone()).expect("Error opening vcf file.");
    let header = bcf.header().clone();
    for line in read_to_string(args.region.clone()).unwrap().lines() {
        regions.push(parse_region(line.to_string(), &header));
    }
    regions
}

fn main() {
    env_logger::init();
    let args = Args::parse();
    let regions = parse_regions(&args);

    let ovl_fn: String = format!("{}.ovls.txt", args.prefix);
    let mut ovl_file = File::create(ovl_fn).expect("Unable to output seq file");

    ovl_file
        .write("#chr\tstart\tref_allele\talt_allele\toverlap_type\tidx\thap\n".as_bytes())
        .unwrap();

    let seq_results_fn: String = format!("{}.haps.fasta", args.prefix);
    let mut seq_file = File::create(seq_results_fn).expect("Unable to output seq file");

    let json_fn: String = format!("{}.variants.json", args.prefix);
    let json_fh = File::create(json_fn).unwrap();
    let mut json_writer = BufWriter::new(json_fh);

    let bcf = &mut IndexedReader::from_path(args.vcf.clone()).expect("Error opening vcf file.");

    let header = bcf.header().clone();
    let sample_count = usize::try_from(header.sample_count()).expect("failure to get samples");

    let mut sample_lookup: HashMap<String, usize> = HashMap::new();
    for sample_index in 0..sample_count {
        let sample_name = String::from_utf8(header.samples()[sample_index].to_vec()).unwrap();
        sample_lookup.insert(sample_name, sample_index);
    }

    let mut json_stuct: HashMap<String, VarTainer> = HashMap::new();

    for r in regions {
        let data = load_data(&args, bcf, &r, &sample_lookup);
        info!("Loaded data for region: {}", r);

        //print!("{:?}\n\n", data.variants[0]);
        //print!("{:?}\n\n", data.variants[1]);

        for h in &data.variants {
            for v in h {
                if v.ovl != OvlType::NoOvl {
                    ovl_file.write(format!("{}\n", v).as_bytes()).unwrap();
                }
            }
        }
        if args.discovery_only == true {
            continue;
        }

        // iterating over both phases
        for (iidx, i) in data.variants.iter().enumerate() {
            info!("Building all unique paths through overlapping variants");
            let all_haps = find_all_paths(&i, 0);
            // iterate over the possible unique haplotype paths
            for (vidx, v) in all_haps.iter().enumerate() {
                let hap = build_haplotypes(&v, &data.sequence, data.region.start as i64);

                let meta = format!("{};{};hap:{}.{}", args.sample.clone(), r, iidx, vidx,);

                let vus: VarTainer = VarTainer {
                    label: meta.clone(),
                    count: v.len(),
                    vars: v.to_vec(),
                };

                json_stuct.insert(meta.clone(), vus);

                seq_file
                    .write(format!(">{} nvar:{}\n{}\n", meta, v.len(), hap).as_bytes())
                    .unwrap();
            }
        }
    }

    json_writer
        .write_all(
            serde_json::to_string_pretty(&json_stuct)
                .unwrap()
                .as_bytes(),
        )
        .unwrap();
    json_writer.write("\n".as_bytes()).unwrap();
    json_writer.flush().unwrap();
    info!("done writing variants to json");
}
