use std::collections::HashMap;
use std::fs::File;
use std::io::Read as IoRead;
use std::path::Path;

use clap::Parser;
use log::LevelFilter;
use rust_htslib::bcf::{Format, Header, Read, Reader, Writer};
use serde::{Deserialize, Serialize};

/// Annotate overlapping variants
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Json (annotation of overlapping variants produced by ovlfilter)
    #[arg(short, long)]
    json: String,

    /// VCF file
    #[arg(short, long)]
    vcf: String,

    /// Output prefix
    #[arg(short, long)]
    prefix: String,

    /// Verbosity
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
}

#[derive(Serialize, Deserialize)]
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

    let mut annotation_info: String = String::new();
    let mut annotation_json = File::open(args.json).unwrap();
    annotation_json
        .read_to_string(&mut annotation_info)
        .unwrap();

    let annotations: Vec<Anno> = serde_json::from_str(&annotation_info).unwrap();
    let mut anno_lookup: HashMap<String, String> = HashMap::new();
    for a in annotations {
        anno_lookup.insert(a.key.clone(), a.anno.clone());
        println!("lu:{} {}", a.key, a.anno);
    }

    let mut bcf = Reader::from_path(args.vcf).expect("Error opening vcf file.");
    let header = bcf.header().clone();
    let mut wheader: Header = Header::from_template(&header);
    wheader.push_record(br#"##INFO=<ID=OVL,Number=1,Type=String,Description="pass/fail info for overlapping variants on the same haplotype">"#);

    let mut outvcf = Writer::from_path(
        Path::new(&format!("{}.ovl_annotated.vcf", args.prefix)),
        &wheader,
        true,
        Format::Vcf,
    )
    .unwrap();
    // 1 overlaps annotation same hap
    // 2 - overlaps annotation same hap
    // json: key allele, value annotation. Two kidns: (overlap [two phased variants overlapping same hap: impossible], overlap but other hap = pass)
    // problem:

    for record_result in bcf.records() {
        let record = record_result.expect("Fail to read record");
        let chr = std::str::from_utf8(header.rid2name(record.rid().unwrap()).expect("Invalid rid"))
            .unwrap();

        let ref_allele: String = (*std::str::from_utf8(record.alleles()[0]).unwrap()).to_string();
        let alt_allele: String = (*std::str::from_utf8(record.alleles()[1]).unwrap()).to_string();

        let lk = format!("{}:{}:{}:{}", chr, record.pos(), ref_allele, alt_allele);

        println!("lk:{}", lk);

        let mut new_record = record;
        outvcf.translate(&mut new_record);
        if anno_lookup.contains_key(&lk) {
            println!("YES");
            let payload = anno_lookup.get(&lk).unwrap().as_bytes();
            new_record.push_info_string(b"OVL", &[payload]).unwrap();
        }
        println!("v: {:?}", new_record);

        outvcf.write(&new_record).unwrap();
    }
}
