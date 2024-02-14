use clap::Parser;
use log::info;
use std::collections::HashSet;
use std::fmt;
use std::fs::read_to_string;
use std::fs::File;
use std::io::Write;

use log::LevelFilter;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    mashed: String,
    #[arg(short, long)]
    samples: String,
    #[arg(short, long)]
    ik: String,
    #[arg(short, long)]
    prefix: String,
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

fn main() {
    let mut failed_sites: i64 = 0;
    let mut passed_sites: i64 = 0;

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

    let mut vcf = File::create(format!("{}.vcf", args.prefix)).expect("Unable to create vcf file");
    let mut fail_bed =
        File::create(format!("{}.fail.bed", args.prefix)).expect("Unable to create bed file");
    let mut pass_bed =
        File::create(format!("{}.pass.bed", args.prefix)).expect("Unable to create bed file");
    let mut upset_r =
        File::create(format!("{}.upsetR.txt", args.prefix)).expect("Unable to create bed file");

    let mut header = "##fileformat=VCFv4.2\n".to_string();
    header.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    header.push_str(
        "##FILTER=<ID=gt_mismatch,Description=\"One or more caller has different genotypes\">\n",
    );
    header.push_str("##INFO=<ID=SOURCES,Number=.,Type=String,Description=\"List of tools or technologies that called the same record\">\n");
    header.push_str("##INFO=<ID=SC,Number=1,Type=Integer,Description=\"Count of supporting tools or technologies\">\n");
    header.push_str(
        "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Allele Type (SNV or INDEL)\">\n",
    );
    header.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    header.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");

    let samples: Vec<String> = args.samples.split(',').map(str::to_string).collect();
    let inputkeys: Vec<String> = args.ik.split(',').map(str::to_string).collect();

    upset_r
        .write(b"site\t")
        .expect("Issue writing to upsetr file");
    for i in &inputkeys {
        upset_r
            .write(format!("{}\t", i).as_bytes())
            .expect("Issue writing to upsetr file");
    }
    upset_r.write(b"\n").expect("Issue writing to upsetr file");

    header.push_str(&samples.join("\t"));
    header.push('\n');
    vcf.write(header.as_bytes()).expect("Failure to write vcf");

    for line in read_to_string(args.mashed).unwrap().lines() {
        let fields: Vec<String> = line.split_whitespace().map(str::to_string).collect();
        let geno: Vec<String> = fields[1].split(',').map(str::to_string).collect();
        let mut sources: Vec<String> = fields
            .last()
            .unwrap()
            .split(',')
            .map(str::to_string)
            .collect();
        sources.sort();

        let source_lookup: HashSet<String> = HashSet::from_iter(sources.iter().cloned());

        upset_r
            .write(format!("{}\t", fields[0]).as_bytes())
            .expect("Issue writing to upsetr file");

        for i in &inputkeys {
            if source_lookup.contains(i) {
                upset_r.write(b"1\t").expect("Issue writing to upsetr file");
            } else {
                upset_r.write(b"0\t").expect("Issue writing to upsetr file");
            }
        }
        upset_r.write(b"\n").unwrap();

        let fg = geno[0].clone();

        let mut seen_genos: HashSet<String> = HashSet::new();

        for g in geno {
            seen_genos.insert(g.clone());
        }

        let mut vcf_filter = "PASS".to_string();

        if seen_genos.len() > 1 {
            failed_sites += 1;
            vcf_filter = "gt_mismatch".to_string();
        }
        passed_sites += 1;

        let pos_allele: Vec<String> = fields[0].split(':').map(str::to_string).collect();

        let mut record_type = "SNV".to_string();
        if pos_allele[3].len() > 1 || pos_allele[4].len() > 1 {
            record_type = "INDEL".to_string();
        }

        let mut record = VcfRecord {
            chrom: pos_allele[0].clone(),
            start: pos_allele[1].clone().parse::<i64>().unwrap() + 1,
            id: "0".to_string(),
            ref_allele: pos_allele[3].clone(),
            alt_allele: pos_allele[4].clone(),
            qual: ".".to_string(),
            filter: vcf_filter.clone(),
            info: format!(
                "TYPE={};SOURCES={};SC={}",
                record_type,
                sources.join(","),
                sources.len()
            ),
            format: "GT".to_string(),
            genotypes: fg.split(':').map(str::to_string).collect(),
        };
        record.genotypes.pop();

        if record.filter == "PASS" {
            pass_bed
                .write(
                    format!("{}\t{}\t{}\n", pos_allele[0], pos_allele[1], pos_allele[2]).as_bytes(),
                )
                .expect("failed to write bed");
        } else {
            fail_bed
                .write(
                    format!("{}\t{}\t{}\n", pos_allele[0], pos_allele[1], pos_allele[2]).as_bytes(),
                )
                .expect("failed to write bed");
        }

        vcf.write(format!("{}\n", record).as_bytes())
            .expect("Failure to write record to vcf");
    }

    info!(
        "passed site count:{} failed site count:{}",
        passed_sites, failed_sites
    );
}
