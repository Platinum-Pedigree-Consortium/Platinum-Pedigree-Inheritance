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
    #[arg(short = 'V', long, help = "Input {b,v}cf path")]
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

/*
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
*/

struct VcfFields {
    // "END=373973;AN=7;NS=7;NA=2;ALEN=74,0;AC=2;VS=>s24;VE=>s25;AWALK=>s67628,*"
    pub end: Option<i32>,
    pub alen: Vec<i32>,
    pub allele_walk: Vec<String>,
    pub vertex_start: String,
    pub vertex_end: String,
}

impl VcfFields {
    pub fn new(record: &bcf::Record) -> Result<Self, Box<dyn std::error::Error>> {
        let end = record.info(b"END"); //("END tag missing");
        let alen = record.info(b"ALEN"); //("ALEN tag missing");
        let awalk = record.info(b"AWALK"); //("AWALK tag missing");
        let vertex_start = record.info(b"VS"); //("VS tag missing");
        let vertex_end = record.info(b"VE"); //("VE tag missing");

        let alen = alen.integer()?.map(|x| x.to_owned()).unwrap_or_default();
        let end = end.integer()?.and_then(|end| end.get(0).copied());
        let awalk = awalk.string()?;
        let allele_walk = awalk
            .map(|awalk| {
                awalk
                    .iter()
                    .map(|awalk| {
                        std::str::from_utf8(&awalk[..])
                            .expect("walk piece was not utf-8")
                            .to_string()
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let vertex_start =
            std::str::from_utf8(vertex_start.string()?.expect("VS was not present")[0])?.to_owned();
        let vertex_end =
            std::str::from_utf8(vertex_end.string()?.expect("VE was not present")[0])?.to_owned();
        Ok(Self {
            end,
            alen,
            allele_walk,
            vertex_start,
            vertex_end,
        })
    }
}

fn core(args: &Args) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("Loading gfa: {:?}.", args.gfa);
    let gfa = gfa::File::from_path(&args.gfa)?;

    /*
    if log::log_enabled!(log::Level::Trace) {
        for (segment_id, segment) in gfa.segments.iter().enumerate() {
            let range = gfa.out_edges(segment_id);
        }
    }
    */

    log::info!("Loaded gfa: {:?}.", args.gfa);
    log::info!("Opening vcf: {:?}.", args.vcf);
    let mut reader = bcf::Reader::from_path(&args.vcf)?;
    let header = bcf::Header::from_template(reader.header()); // TODO: add new header info potentially. Maybe auto-handle missing vcf fields.

    let bcf = args
        .bcf
        .as_ref()
        .map_or("-".to_string(), |x| x.display().to_string());
    let mut writer = bcf::Writer::from_path(&bcf, &header, args.uncompressed, bcf::Format::Bcf)?;
    log::info!("Opened handles; reading in data.");
    let mut all_records = Vec::new();
    //let mut reader = bcf::Reader::from_path(&args.vcf)?;
    for record in reader.records() {
        let mut record = record?;
        writer.translate(&mut record);
        record.unpack();
        let record_data = VcfFields::new(&record)?;
        all_records.push((record, record_data));
    }
    log::info!("Read in data, processing");
    for (record, data) in &mut all_records {
        let vertex_start = gfa.segment_id(&data.vertex_start).expect("No VS");
        let vertex_end = gfa.segment_id(&data.vertex_end).expect("No VE");
        log::debug!("VE,VS: {vertex_start}, {vertex_end}");
    }

    log::info!("Processed, writing out");
    for (record, _data) in &all_records {
        writer.write(record)?;
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
