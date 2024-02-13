use clap::Parser;
use concordance::gfa;
use log::info;
use std::boxed::Box;
use std::collections::HashSet;
use std::fmt;
use std::fs::read_to_string;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use log::LevelFilter;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[clap(help = "Path to input gfa")]
    gfa: PathBuf,
    #[arg(short, long, help = "Input vcf path")]
    vcf: PathBuf,
    #[arg(short, long, help = "Output bcf path")]
    bcf: PathBuf,
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
}

#[derive(Debug)]
struct VcfRecord {
    chrom: String,
    start: i64,
    id: String,
    ref_allele: String,
    alt_allele: String,
    qual: String,
    filter: String,
    info: String,
    format: String,
    genotypes: Vec<String>,
}

impl fmt::Display for VcfRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.chrom,
            self.start,
            self.id,
            self.ref_allele,
            self.alt_allele,
            self.qual,
            self.filter,
            self.info,
            self.format,
            self.genotypes.join("\t")
        )
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
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

    let gfa = gfa::File::from_path(&args.gfa)?;
    if log::log_enabled!(log::Level::Debug) {
        for (segment_id, segment) in gfa.segments.iter().enumerate() {
            let range = gfa.out_edges(segment_id);
            eprintln!("segment {} has {} out edges", segment.name(), range.count());
        }
    }
    log::info!("Loaded gfa: {gfa}");
    Ok(())
}
