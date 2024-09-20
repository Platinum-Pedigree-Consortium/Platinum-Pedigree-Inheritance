use std::fs::File;
use std::io::Read as IoRead;
use std::path::Path;
use std::{collections::HashMap, fmt::Debug};

use clap::Parser;
use log::LevelFilter;
use rust_htslib::bcf::{Format, Header, Read, Reader, Writer};
use serde::{Deserialize, Serialize};
use serde_json::to_string;

/// Annotate overlapping variants
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Json (annotation of overlapping variants produced by ovlfilter)
    #[arg(short, long)]
    json: String,

    /// VCF file
    #[arg(short, long)]
    data: String,

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
    let mut annotation_json = File::open(&args.json)
        .unwrap_or_else(|e| panic!("Failed to open file at {:?}. Error: {e:?}", args.json));
    annotation_json
        .read_to_string(&mut annotation_info)
        .unwrap();

    let annotations: Vec<Anno> = serde_json::from_str(&annotation_info).unwrap();
    let mut anno_lookup: HashMap<String, String> = HashMap::new();
    for a in annotations {
        anno_lookup.insert(a.key.clone(), a.anno.clone());
    }

    let mut bcf = Reader::from_path(args.data).expect("Error opening vcf file.");
    let header = bcf.header().clone();
    let mut wheader: Header = Header::from_template(&header);
    wheader.remove_info(b"OVL");
    wheader.push_record(br#"##INFO=<ID=OVL,Number=.,Type=String,Description="annotation for overlapping variants on the same haplotype">"#);

    let mut outvcf = Writer::from_path(
        Path::new(&format!("{}.ovl_annotated.vcf", args.prefix)),
        &wheader,
        true,
        Format::Vcf,
    )
    .unwrap();

    for (record_number, record_result) in bcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let chr = std::str::from_utf8(header.rid2name(record.rid().unwrap()).expect("Invalid rid"))
            .unwrap();

        let ref_allele: String = (*std::str::from_utf8(record.alleles()[0]).unwrap()).to_string();
        let alt_allele: String = (*std::str::from_utf8(record.alleles()[1]).unwrap()).to_string();

        let lk = format!("{}:{}:{}:{}", chr, record.pos(), ref_allele, alt_allele);

        let mut new_record = record;
        let name = String::from_utf8(new_record.id()).unwrap();

        outvcf.translate(&mut new_record);
        if let Some(payload) = anno_lookup.get(&lk).map(|x| x.as_bytes()) {
            let mut source_string = if let Ok(Some(string)) = new_record.info(b"OVL").string() {
                let slice = &string[..];
                let seqs = slice
                    .iter()
                    .map(|x| std::str::from_utf8(x).unwrap())
                    .map(std::borrow::ToOwned::to_owned)
                    .collect::<Vec<_>>();
                seqs
            } else {
                vec![]
            };
            log::trace!(
                "OVL string before: {source_string:?} for record {name} at count {record_number}",
            );
            source_string.push(String::from_utf8(payload.to_owned()).unwrap());
            log::trace!(
                "OVL string after adding payload: {:?}",
                source_string
                    .iter()
                    .map(|x| format!("{x:?}")) // std::str::from_utf8(x).map(std::borrow::ToOwned::to_owned).unwrap_or_else(|_| format!("{x:?}")))
                    .collect::<Vec<_>>()
            );
            let new_tag = source_string
                .iter()
                .map(|x| x.as_bytes())
                .collect::<Vec<_>>();
            new_record.push_info_string(b"OVL", &new_tag[..]).unwrap();
        }
        outvcf.write(&new_record).unwrap();
    }
}
