use concordance::gfa;

use clap::Parser;
use log::LevelFilter;
use rust_htslib::bcf::{self, Read};

use std::boxed::Box;
use std::fmt;
use std::path::PathBuf;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[clap(help = "Path to input gfa")]
    gfa: PathBuf,
    #[arg(short, long, help = "Input {b,v}cf path")]
    vcf: PathBuf,
    #[arg(short, long, help = "Output bcf path")]
    bcf: Option<PathBuf>,
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
    #[arg(
        short,
        long,
        help = "Emit uncompressed output",
        default_value_t = false
    )]
    uncompressed: bool,
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

fn core(args: &Args) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("Loading gfa: {:?}.", args.gfa);
    let gfa = gfa::File::from_path(&args.gfa)?;

    if log::log_enabled!(log::Level::Debug) {
        for (segment_id, segment) in gfa.segments.iter().enumerate() {
            let range = gfa.out_edges(segment_id);
            eprintln!("segment {} has {} out edges", segment.name(), range.count());
        }
    }

    log::info!("Loaded gfa: {:?}.", args.gfa);
    log::info!("Opening vcf: {:?}.", args.vcf);
    let mut reader = bcf::Reader::from_path(&args.vcf)?;
    let header = bcf::Header::from_template(reader.header()); // TODO: add new header info potentially.

    let bcf = args
        .bcf
        .as_ref()
        .map_or("-".to_string(), |x| x.display().to_string());
    let mut writer = bcf::Writer::from_path(&bcf, &header, args.uncompressed, bcf::Format::Bcf)?;
    for record in reader.records() {
        let mut record = record?;
        writer.translate(&mut record);
        // Process
        writer.write(&record)?;
    }

    Ok(())
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

    core(&args)?;

    Ok(())
}
