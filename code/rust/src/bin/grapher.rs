use std::collections::HashMap;

use std::{fmt, fs::File, io::BufWriter, io::Write, str};

use clap::Parser;

use rust_htslib::bcf::{record::GenotypeAllele, IndexedReader, Read};

use concordance::vcf2seq::build_haplotypes;
use concordance::vcf2seq::find_all_paths;
use concordance::vcf2seq::VarTainer;

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
