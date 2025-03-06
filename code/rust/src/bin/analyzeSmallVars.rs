use clap::Parser;

use rust_htslib::bcf::{Read, Reader};
use std::collections::HashSet;

/// A tool to extract info for analyzing/plotting tool overlaps
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// VCF file
    #[arg(short, long)]
    vcf: String,
    /// list of "SOURCES" in info field to count
    #[arg(short, long)]
    supports: String,
}

fn main() {
    let args = Args::parse();
    let mut bcf = Reader::from_path(args.vcf).expect("Error opening vcf file.");

    let tools: Vec<String> = args.supports.split(",").map(str::to_string).collect();

    println!("{}", tools.join("\t"));

    for (_i, record_result) in bcf.records().enumerate() {
        let mut sources: HashSet<String> = HashSet::new();
        let record = record_result.expect("Fail to read record");
        let b_sources = record.info(b"SOURCES").string().unwrap().unwrap();

        for s in b_sources.iter() {
            sources.insert((*std::str::from_utf8(*s).unwrap()).to_string());
        }

        let mut rez: Vec<String> = Vec::new();

        for t in &tools {
            if sources.contains(t) {
                rez.push("1".to_string());
            } else {
                rez.push("0".to_string());
            }
        }

        println!("{}", rez.join("\t"));
    }
}
