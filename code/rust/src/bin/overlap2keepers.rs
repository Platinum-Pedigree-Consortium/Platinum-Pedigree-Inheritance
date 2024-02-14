use std::collections::{HashMap, HashSet};
use std::fs::read_to_string;
use std::i32;
use std::io::Read;
use std::{fs::File, io::Write, str};

use clap::Parser;
use log::{info, warn, LevelFilter};
use serde::{Deserialize, Serialize};

/// A tool for sorting which haplotypes to keep
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// minimap2 paf with alignment
    #[arg(short, long)]
    paf: String,

    /// variant json, mapping from variants to haplotype
    #[arg(short, long)]
    variants: String,

    /// prefix for output
    #[arg(short, long)]
    prefix: String,

    // min alignment block length (paf col. 11)
    #[arg(short, long, default_value_t = 1500)]
    block: i32,

    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 1)]
    verbosity: u8,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
enum VarType {
    Snv,
    Insertion,
    Deletion,
    Ref,
}

#[derive(Debug, PartialEq, Eq, Clone, Serialize, Deserialize)]
enum OvlType {
    NoOvl,
    OvlNext,
    OvlPrev,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
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
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct VarTainer {
    label: String,
    count: usize,
    vars: Vec<Var>,
}

#[derive(Serialize)]
struct Anno {
    key: String,
    anno: String,
}

fn main() {
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

    let mut lowest_edits: HashMap<String, u16> = HashMap::new();
    let mut lowest_ed_hap: HashMap<String, String> = HashMap::new();

    let mut variant_info: String = String::new();

    let mut variant_json = File::open(args.variants).unwrap();
    let mut varriant_anno_fh = File::create(format!("{}.var_anno.json", args.prefix)).unwrap();
    let mut variant_anno: Vec<Anno> = Vec::new();
    let mut variant_marker: HashSet<String> = HashSet::new();

    variant_json.read_to_string(&mut variant_info).unwrap();

    let haplotype_vars: HashMap<String, VarTainer> = serde_json::from_str(&variant_info).unwrap();
    info!("num vars {} ", haplotype_vars.len());

    for line in read_to_string(args.paf.clone()).unwrap().lines() {
        let paf_parts: Vec<String> = line.split('\t').map(str::to_string).collect();
        let edit_dist_parts: Vec<String> = paf_parts
            .get(12)
            .unwrap()
            .split(':')
            .map(str::to_string)
            .collect();
        let hapkey = paf_parts.first().unwrap().clone();
        let hapkey_parts: Vec<String> = hapkey.split('.').map(str::to_string).collect();
        let final_key = hapkey_parts.first().unwrap();
        let aln_len = paf_parts.get(10).unwrap().parse::<i32>().unwrap();

        let edit_distance = edit_dist_parts.get(2).unwrap().parse::<u16>().unwrap();
        lowest_edits
            .entry(final_key.to_string())
            .or_insert(edit_distance);

        lowest_ed_hap
            .entry(final_key.to_string())
            .or_insert(hapkey.clone());

        let v = lowest_edits.get(&final_key.to_string()).unwrap();

        if aln_len < args.block {
            warn!("skiping entry with low block len: {}", final_key);
            continue;
        }

        if edit_distance < *v {
            info!("prev: {} better: {} key: {}", v, edit_distance, final_key);
            *lowest_edits.get_mut(&final_key.to_string()).unwrap() = edit_distance;
            *lowest_ed_hap.get_mut(&final_key.to_string()).unwrap() = hapkey.clone();
        }
    }

    for (i, _j) in lowest_edits {
        let k = lowest_ed_hap.get(&i).unwrap();
        for v in &haplotype_vars.get(k).unwrap().vars {
            let vk = format!("{}:{}:{}:{}", v.seqid, v.start, v.ref_allele, v.alt_allele);
            variant_anno.push(Anno {
                key: vk.clone(),
                anno: "HAP_SELECTED".to_string(),
            });
            variant_marker.insert(vk);
        }
    }

    // going over all the variants to find variants that should be tossed
    for (_k, v) in haplotype_vars {
        for v in v.vars {
            let vk = format!("{}:{}:{}:{}", v.seqid, v.start, v.ref_allele, v.alt_allele);
            if !variant_marker.contains(&vk) {
                variant_anno.push(Anno {
                    key: vk,
                    anno: "HAP_DESELECTED".to_string(),
                })
            }
        }
    }

    varriant_anno_fh
        .write(
            serde_json::to_string_pretty(&variant_anno)
                .unwrap()
                .as_bytes(),
        )
        .expect("Failure to write variant json");
    varriant_anno_fh
        .write("\n".as_bytes())
        .expect("Failure to write to variant json");
}
