use concordance::gfa::{self, Orientation};

use clap::Parser;
use itertools::Itertools;
use log::LevelFilter;
use rust_htslib::bcf::{self, Read};

use std::boxed::Box;
use std::path::PathBuf;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[clap(help = "Path to input gfa")]
    gfa: PathBuf,
    #[arg(help = "Input {b,v}cf path")]
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

pub fn rc(x: &str) -> String {
    x.chars()
        .rev()
        .map(|x| match x {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            't' => 'a',
            _ => '?',
        })
        .collect::<String>()
}

#[derive(Debug)]
struct VariantDatum {
    // "END=373973;AN=7;NS=7;NA=2;ALEN=74,0;AC=2;VS=>s24;VE=>s25;AWALK=>s67628,*"
    pub alen: Vec<i32>,
    pub allele_walk: Vec<String>,
    vertex_start: String,
    vertex_end: String,
}

impl VariantDatum {
    pub fn new(record: &bcf::Record) -> Result<Self, Box<dyn std::error::Error>> {
        //let end = record.info(b"END"); //("END tag missing");
        let alen = record.info(b"ALEN"); //("ALEN tag missing");
        let awalk = record.info(b"AWALK"); //("AWALK tag missing");
        let vertex_start = record.info(b"VS"); //("VS tag missing");
        let vertex_end = record.info(b"VE"); //("VE tag missing");

        let alen = alen.integer()?.map(|x| x.to_owned()).unwrap_or_default();
        //let end = end.integer()?.and_then(|end| end.first().copied());
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
            //end,
            alen,
            allele_walk,
            vertex_start,
            vertex_end,
        })
    }

    /// Get id for vertex for start.
    /// Trims the leading '>' character.
    /// If no '>' is present, uses the string directly.
    pub fn vertex_start(&self) -> &str {
        let start = usize::from(self.vertex_start.starts_with('>'));
        &self.vertex_start[start..]
    }

    /// Get id for vertex for end.
    /// Trims the leading '>' character.
    /// If no '>' is present, uses the string directly.
    pub fn vertex_end(&self) -> &str {
        let start = usize::from(self.vertex_end.starts_with('>'));
        &self.vertex_end[start..]
    }
}

/*
#[derive(Debug)]
struct WalkStep(pub String, pub Orientation);
*/

fn extract_seq(input: (&Orientation, &String), gfa: &gfa::File) -> String {
    let (orient, vertex) = input;
    if vertex == "*" {
        return "".into();
    }
    let vtx = gfa
        .segment_id(&vertex)
        .unwrap_or_else(|| panic!("Missing vertex {vertex}"));
    let vtx = &gfa.segments[vtx];
    let seq = vtx
        .sequence
        .as_ref()
        .expect("segment for vertex does not have a sequence");
    if *orient == Orientation::Forward {
        seq.clone()
    } else {
        rc(&seq)
    }
}

fn dot_if_empty(mut x: String) -> String {
    if x.is_empty() {
        x.push('.');
    }
    x
}

fn make_seq(alleles: &[(Orientation, String)], gfa: &gfa::File) -> String {
    dot_if_empty(
        Itertools::intersperse(
            alleles.iter().map(|(orient, seq)| match orient {
                Orientation::Forward | Orientation::Reverse => extract_seq((orient, seq), gfa),
                Orientation::Star => "".into(),
            }),
            "".to_string(),
        )
        .collect::<String>(),
    )
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
    let reader_for_header = bcf::Reader::from_path(&args.vcf)?;
    let header_view = reader_for_header.header();

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
        let record_data = VariantDatum::new(&record)?;
        all_records.push((record, record_data));
    }
    log::info!("Read in data, processing");
    // chr1:148288481
    // chr1:150411409
    const _TRUE_INSERTION: (i32, i32) = (0, 148_288_482 - 1); //chr1    148288482
    const _TRUE_DELETION: (i32, i32) = (0, 150_411_410 - 1); // chr1    150411410
    const _TRUTH: &[i32] = &[_TRUE_INSERTION.1, _TRUE_DELETION.1];
    for (record, data) in &mut all_records {
        let vertex_start_str = data.vertex_start();
        log::trace!("Getting vertex start: {vertex_start_str}");
        let vertex_start = gfa.segment_id(&vertex_start_str).unwrap_or_else(|| {
            panic!(
                "No vertex found for data {data:?} with vs {vertex_start_str}/{:?}",
                data.vertex_start
            )
        });
        let vertex_end_str = data.vertex_end();
        log::trace!("Getting vertex end: {vertex_end_str}");
        let vertex_end = gfa.segment_id(&vertex_end_str).unwrap_or_else(|| {
            panic!(
                "No vertex found for data {data:?} with ve {:?}",
                data.vertex_end
            )
        });
        log::trace!("VE,VS: {vertex_start}, {vertex_end} ({vertex_start_str}, {vertex_end_str})");
        let walks = data
            .allele_walk
            .iter()
            .map(|x| {
                if x == "*" {
                    vec![(Orientation::Star, "*".to_owned())]
                } else {
                    gfa::decompose_walk(x)
                }
            })
            .collect::<Vec<Vec<(Orientation, String)>>>();
        assert_eq!(walks.len(), data.alen.len());
        let seqs = walks
            .iter()
            .map(|x| make_seq(&x[..], &gfa))
            .collect::<Vec<_>>();
        let alens = seqs
            .iter()
            .map(|x| if x == "." { 0 } else { x.len() as i32 })
            .collect::<Vec<_>>();
        assert_eq!(
            alens, data.alen,
            "seqs: {seqs:?}. expected alens: {:?}, found {alens:?}",
            data.alen
        );
        let id = String::from_utf8(record.id()).expect("id not utf-8");
        let id = if id == "." {
            let length_str =
                Itertools::intersperse(alens.iter().map(|x| x.to_string()), ",".to_string())
                    .collect::<String>();
            format!(
                "{}:{}:{length_str}",
                record
                    .rid()
                    .and_then(|x| header_view.rid2name(x).ok())
                    .map_or(record.rid().map_or("*".into(), |x| x.to_string()), |x| {
                        String::from_utf8(x.to_owned())
                            .expect("Failed to utf-8 encode a reference name")
                    }),
                record.pos()
            )
        } else {
            id
        };
        record.set_id(id.as_bytes())?;
        if !seqs.is_empty() {
            let alleles = seqs
                .iter()
                .map(|x| &x.as_bytes()[..])
                .collect::<Vec<&[u8]>>();
            record.set_alleles(&alleles[..])?;
        }
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
