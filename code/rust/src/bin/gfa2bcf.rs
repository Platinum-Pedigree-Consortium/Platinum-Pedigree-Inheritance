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
REF - reference base(s): Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted.
The value in the POS field refers to the position of the first base in the String.
---
For simple insertions and
deletions in which either the REF or one of the ALT alleles would otherwise be null/empty, the REF and ALT
Strings must include the base before the event (which must be reflected in the POS field), unless the event
occurs at position 1 on the contig in which case it must include the base after the event; this padding base is
not required (although it is permitted) for e.g. complex substitutions or other events where all alleles have at
least one base represented in their Strings.
---
If any of the ALT alleles is a symbolic allele (an angle-bracketed
ID String “<ID>”) then the padding base is required and POS denotes the coordinate of the base preceding
the polymorphism. Tools processing VCF files are not required to preserve case in the allele Strings. (String,
Required).
5. ALT - alternate base(s): Comma separated list of alternate non-reference alleles. These alleles do not have to
be called in any of the samples. Options are base Strings made up of the bases A,C,G,T,N,*, (case insensitive)
or an angle-bracketed ID String (“<ID>”) or a breakend replacement string as described in the section on
breakends. The ‘*’ allele is reserved to indicate that the allele is missing due to a upstream deletion. If there
are no alternative alleles, then the missing value should be used. Tools processing VCF files are not required
to preserve case in the allele String, except for IDs, which are case sensitive. (String; no whitespace, commas,
or angle-brackets are permitted in the ID String itself)
*/

#[inline]
pub fn rc_char(x: char) -> char {
    match x {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        'a' => 't',
        'c' => 'g',
        'g' => 'c',
        't' => 'a',
        _ => '?',
    }
}

pub fn rc(x: &str) -> String {
    x.chars().rev().map(rc_char).collect::<String>()
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

    pub fn vertex_start_orientation(&self) -> Orientation {
        Orientation::new(self.vertex_start.chars().next().unwrap())
    }
    /*
    pub fn vertex_end_orientation(&self) -> Orientation {
        Orientation::new(self.vertex_end.chars().next().unwrap())
    }
    */

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

fn star_if_empty(mut x: String) -> String {
    if x.is_empty() {
        x.push('*');
    }
    x
}

fn make_seq(alleles: &[(Orientation, String)], gfa: &rgfa::File) -> String {
    star_if_empty(
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
    {
        if walks.iter().map(|x| &x[0].1).unique().count() == 1 {
            assert_eq!(
            walks.iter().map(|x| x[0].0).unique().count(),
            2,
            "Make sure that the two identical walks are used in opposite orientations. Walks: {walks:?}"
        );
            // If the same walk is used in both alleles, then
            "INV"
        } else {
            "MNP" // Not an inversion, same length.
        }
    } else {
        "SV"
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
        .map_or(false, |x| match &x[..] {
            "vcf" => false,
            "bcf" => true,
            _ => {
                log::warn!("Unexpected extension {x:?}. Defaulting to vcf format.)");
                false
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
        let initial_pos = record.pos();
        assert_eq!(walks.len(), data.alen.len());
        let mut seqs = walks
            .iter()
            .map(|x| make_seq(&x[..], &gfa))
            .collect::<Vec<_>>();
        let alens = seqs
            .iter()
            .map(|x| if x == "*" { 0 } else { x.len() as i32 })
            .collect::<Vec<_>>();
        let reflen_zero = alens.first() == Some(&0);
        let varlen_zero = alens[1..].iter().min() == Some(&0);
        // Pad with previous base if necessary.
        /*
        VCF4.2 spec: if any allele is empty (IE, min length is 0), then these strings must include the base before the event.
        ```
        For simple insertions and
        deletions in which either the REF or one of the ALT alleles would otherwise be null/empty, the REF and ALT
        Strings must include the base before the event (which must be reflected in the POS field), unless the event
        occurs at position 1 on the contig in which case it must include the base after the event; this padding base is
        not required (although it is permitted) for e.g. complex substitutions or other events where all alleles have at
        least one base represented in their Strings.
        ```
        */
        let either_zero = reflen_zero || varlen_zero;
        // We have to subtract one from position if we are adding a reference base here.
        let final_pos = initial_pos - i64::from(either_zero);
        if either_zero {
            let preceding_seq = gfa.segments[vertex_start]
                .sequence
                .as_ref()
                .ok_or_else(|| {
                    rgfa::Error(format!(
                        "segment {vertex_start_str}/{vertex_start} had no sequence"
                    ))
                })?;
            let vertex_start_orientation = data.vertex_start_orientation();
            let padding_char = match vertex_start_orientation {
                Orientation::Forward => preceding_seq.chars().last().unwrap_or('*'),
                Orientation::Reverse => rc_char(preceding_seq.chars().next().unwrap_or('*')),
                Orientation::Star => {
                    return Err(Box::new(rgfa::Error("Unexpected * orientation".into())));
                }
            };
            for seq in &mut seqs {
                *seq = if seq == "*" {
                    padding_char.to_string()
                } else {
                    format!("{padding_char}{seq}")
                }
            }
        }

        // "*" allele if none present.
        if seqs.len() == 1 {
            seqs.push("*".into());
        }

        log::debug!("{} total seqs, seqs {seqs:?}", seqs.len());
        assert_eq!(
            alens, data.alen,
            "seqs: {seqs:?}. expected alens: {:?}, found {alens:?}",
            data.alen
        );
        let mut id = String::from_utf8(record.id()).expect("id not utf-8");
        if id == "." {
            let length_str =
                Itertools::intersperse(alens.iter().map(|x| x.to_string()), ",".to_string())
                    .collect::<String>();
            id = length_str;
        }
        let id = format!("{}:{id}:{final_pos}", id_from_alens(&alens, &walks));
        record.set_id(id.as_bytes())?;
        if final_pos != initial_pos {
            record.set_pos(final_pos);
        }

        // Assign alleles
        if !seqs.is_empty() {
            let alleles = seqs.iter().map(|x| x.as_bytes()).collect::<Vec<&[u8]>>();
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
