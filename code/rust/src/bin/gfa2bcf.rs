use concordance::rgfa::{self, Orientation};

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
    #[arg(short, long, help = "Output [bv]cf path")]
    output: Option<PathBuf>,
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
    #[arg(
        short,
        long,
        help = "Emit uncompressed output. If unset, emits compressed for bcf and uncompressed for vcf."
    )]
    uncompressed: Option<bool>,
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

fn extract_seq(input: (&Orientation, &String), gfa: &rgfa::File) -> String {
    let (orient, vertex) = input;
    if vertex == "*" {
        return "".into();
    }
    let vtx = gfa
        .segment_id(vertex)
        .unwrap_or_else(|| panic!("Missing vertex {vertex}"));
    let vtx = &gfa.segments[vtx];
    let seq = vtx
        .sequence
        .as_ref()
        .expect("segment for vertex does not have a sequence");
    if *orient == Orientation::Forward {
        seq.clone()
    } else {
        rc(seq)
    }
}

fn dot_if_empty(mut x: String) -> String {
    if x.is_empty() {
        x.push('.');
    }
    x
}

fn make_seq(alleles: &[(Orientation, String)], gfa: &rgfa::File) -> String {
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

fn id_from_alens(alens: &[i32], walks: &[Vec<(Orientation, String)>]) -> &'static str {
    if alens.len() <= 1 {
        "."
    } else if alens[0] == 0 && alens[1..].iter().all(|x| *x > 0) {
        "INS"
    } else if alens[0] > 0 && alens[1..].iter().all(|x| *x == 0) {
        "DEL"
    } else if alens.len() == 2
        && alens[0] == alens[1]
        && *alens.iter().max().unwrap() > 0
        && walks.iter().all(|x| x.len() == 1)
        && walks.iter().map(|x| &x[0].1).unique().count() == 1
    {
        assert_eq!(
            walks.iter().map(|x| x[0].0).unique().count(),
            2,
            "Make sure that the two identical walks are used in opposite orientations. Walks: {walks:?}"
        );
        // If the same walk is used in both alleles, then
        "INV"
    } else {
        "."
    }
}

fn core(args: &Args) -> Result<(), Box<dyn std::error::Error>> {
    log::info!("Loading gfa: {:?}.", args.gfa);
    let gfa = rgfa::File::from_path(&args.gfa)?;

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

    let use_bcf = args
        .output
        .as_ref()
        .and_then(|bcf| {
            bcf.extension()
                .map(|x| x.to_str().expect("Could not convert to utf-8").to_string())
        })
        .map_or(true, |x| match &x[..] {
            "vcf" => false,
            "bcf" => true,
            _ => {
                log::warn!("Unexpected extension {x:?}. Defaulting to bcf format.)");
                true
            }
        });
    let format = if use_bcf {
        bcf::Format::Bcf
    } else {
        bcf::Format::Vcf
    };
    let bcf = args
        .output
        .as_ref()
        .map_or("-".to_string(), |x| x.display().to_string());
    let uncompressed = args.uncompressed.unwrap_or(!use_bcf);
    // Uncompressed vcf, compressed bcf by default. Can override with `-u`.

    let mut writer = bcf::Writer::from_path(bcf, &header, uncompressed, format)?;
    log::info!("Opened handles; reading in data.");
    // TODO: rewrite as streaming.
    let mut all_records = Vec::new();
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
        let vertex_start = gfa.segment_id(vertex_start_str).unwrap_or_else(|| {
            panic!(
                "No vertex found for data {data:?} with vs {vertex_start_str}/{:?}",
                data.vertex_start
            )
        });
        let vertex_end_str = data.vertex_end();
        log::trace!("Getting vertex end: {vertex_end_str}");
        let vertex_end = gfa.segment_id(vertex_end_str).unwrap_or_else(|| {
            panic!(
                "No vertex found for data {data:?} with ve {:?}",
                data.vertex_end
            )
        });
        log::trace!("VE,VS: {vertex_start}, {vertex_end} ({vertex_start_str}, {vertex_end_str})");
        let mut walks = Vec::new();
        for walk in &data.allele_walk {
            walks.push(if walk == "*" {
                vec![(Orientation::Star, "*".to_owned())]
            } else {
                rgfa::decompose_walk(walk)?
            })
        }
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
        let pos = record.pos();
        let mut id = String::from_utf8(record.id()).expect("id not utf-8");
        if id == "." {
            let length_str =
                Itertools::intersperse(alens.iter().map(|x| x.to_string()), ",".to_string())
                    .collect::<String>();
            id = length_str;
        }
        let orig_genotypes = format!("{:?}", record.genotypes()?);
        let orig_genotypes = orig_genotypes.split_terminator("Buffer").next().unwrap();
        let id = format!("{}:{id}:{pos}", id_from_alens(&alens, &walks));
        record.set_id(id.as_bytes())?;
        if !seqs.is_empty() {
            let mut alleles: Vec<&[u8]> = vec![b"N"];
            for seq in &seqs {
                alleles.push(seq.as_bytes());
            }
            record.set_alleles(&alleles[..])?;
        }
        let final_genotypes = format!("{:?}", record.genotypes()?);
        let final_genotypes = final_genotypes.split_terminator("Buffer").next().unwrap();
        assert_eq!(
            final_genotypes, orig_genotypes,
            "Old geno: {orig_genotypes}. Now geno: {final_genotypes}"
        );
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
