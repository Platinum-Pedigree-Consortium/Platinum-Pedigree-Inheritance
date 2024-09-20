use clap::Parser;
use log::info;
use rust_htslib::bcf::record::Genotype;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt;
use std::fs::read_to_string;
use std::fs::File;
use std::io::Write;
use std::vec;

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
    #[arg(short, long)]
    fai: String,
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

struct MashRecord {
    chrom: String,
    start: String,
    end: String,
    ref_allele: String,
    alt_alleles: Vec<String>,
    genotype_sets: HashSet<String>,
    sources_lookup: HashSet<String>,
}

impl MashRecord {
    fn genotype_conflict(&self) -> bool {
        if self.genotype_sets.len() > 1 {
            return true;
        }
        false
    }
    fn chr(&self) -> String {
        self.chrom.clone()
    }
    fn pos(&self) -> i64 {
        self.start.clone().parse::<i64>().unwrap() + 1
    }
    fn end(&self) -> i64 {
        self.end.clone().parse::<i64>().unwrap() + 1
    }
    fn genotypes(&self) -> Vec<String> {
        let mut geno_str: Vec<String> = self
            .genotype_sets
            .iter()
            .next()
            .unwrap()
            .split(":")
            .map(str::to_string)
            .collect();
        geno_str.pop();
        geno_str
    }

    fn reference(&self) -> String {
        self.ref_allele.clone()
    }

    fn alts(&self) -> Vec<String> {
        self.alt_alleles.clone()
    }

    fn mash2vcf(&self) -> VcfRecord {
        let sources_vec: Vec<String> = self.sources_lookup.clone().into_iter().collect();

        let mut vrecord = VcfRecord {
            chrom: self.chr(),
            start: self.pos(),
            id: "0".to_string(),
            ref_allele: self.ref_allele.clone(),
            alt_allele: self.alts().join(","),
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            info: format!(
                "SOURCES={};SC={}",
                sources_vec.join(","),
                self.sources_lookup.len()
            ),
            format: "GT".to_string(),
            genotypes: self.genotypes(),
        };
        if self.genotype_conflict() {
            vrecord.filter = "gt_mismatch".to_string();
        }

        vrecord
    }
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

struct Contig {
    id: String,
    len: String,
}

fn parse_contigs(fai: String) -> Vec<Contig> {
    let mut result: Vec<Contig> = Vec::new();
    for line in read_to_string(fai).unwrap().lines() {
        let fields: Vec<String> = line.split_whitespace().map(str::to_string).collect();
        result.push(Contig {
            id: fields.get(0).unwrap().clone(),
            len: fields.get(1).unwrap().clone(),
        });
    }
    return result;
}

fn load_multi(mashed: String) -> HashMap<String, Vec<(i64, String, String)>> {
    let mut result: HashMap<String, Vec<(i64, String, String)>> = HashMap::new();

    let mut line_n = 0;
    for line in read_to_string(mashed).unwrap().lines() {
        line_n += 1;
        let record = parse_record(line);
        if record.genotype_conflict() {
            continue;
        }
        for a in record.alts() {
            let lookup = format!("{}:{}", record.chr(), record.pos());
            result
                .entry(lookup)
                .or_insert(Vec::new())
                .push((line_n, record.reference(), a));
        }
    }
    result
}

fn parse_record(input: &str) -> MashRecord {
    let fields: Vec<String> = input.split_whitespace().map(str::to_string).collect();
    let pos_allele: Vec<String> = fields[0].split(":").map(str::to_string).collect();
    let alts: Vec<String> = pos_allele[4]
        .clone()
        .split(",")
        .map(str::to_string)
        .collect();
    let geno: Vec<String> = fields[1].split(",").map(str::to_string).collect();
    let mut sources: Vec<String> = fields
        .last()
        .unwrap()
        .split(",")
        .map(str::to_string)
        .collect();
    sources.sort();

    let mut results = MashRecord {
        chrom: pos_allele[0].clone(),
        start: pos_allele[1].clone(),
        end: pos_allele[2].clone(),
        ref_allele: pos_allele[3].clone(),
        alt_alleles: alts,
        genotype_sets: HashSet::from_iter(geno.iter().cloned()),
        sources_lookup: HashSet::from_iter(sources.iter().cloned()),
    };
    results
}

fn main() {
    let mut failed_sites: i64 = 0;
    let mut contained_sites: i64 = 0;
    let mut gt_missmatch: i64 = 0;
    let mut passed_sites: i64 = 0;

    let args = Args::parse();

    let filter_level: LevelFilter = match args.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };

    let contigs = parse_contigs(args.fai);

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
    header.push_str(
        "##FILTER=<ID=contained,Description=\"One or more caller has overlapping alleles\">\n",
    );
    header.push_str("##INFO=<ID=SOURCES,Number=.,Type=String,Description=\"List of tools or technologies that called the same record\">\n");
    header.push_str("##INFO=<ID=SC,Number=1,Type=Integer,Description=\"Count of supporting tools or technologies\">\n");
    header.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    for c in contigs {
        header.push_str(&format!("##contig=<ID={},length={}>\n", c.id, c.len));
    }

    header.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");

    let samples: Vec<String> = args.samples.split(",").map(str::to_string).collect();
    let inputkeys: Vec<String> = args.ik.split(",").map(str::to_string).collect();

    header.push_str(&samples.join("\t"));
    header.push_str("\n");
    vcf.write(header.as_bytes()).expect("Failure to write vcf");

    let lookup_alleles = load_multi(args.mashed.clone());

    let mut line_n: i64 = 0;

    for line in read_to_string(args.mashed.clone()).unwrap().lines() {
        line_n += 1;

        let record = parse_record(line);

        let mut vrecord = record.mash2vcf();

        let mut hit_count = 0;

        for a in record.alts() {
            let lookup = format!("{}:{}", record.chr(), record.pos());

            match lookup_alleles.get(&lookup) {
                Some(value) => {
                    for i in value {
                        if value.len() > 1 {}
                        // duplicate entry (by line number), not same reference
                        if i.0 == line_n || record.reference() != i.1 {
                            continue;
                        }
                        // we found a match alt allele
                        if i.2 == a {
                            hit_count += 1;
                        }
                    }
                }
                None => {}
            }
        }

        if record.genotype_conflict() {
            failed_sites += 1;
            gt_missmatch += 1;
            vrecord.filter = "gt_mismatch".to_string();
        }
        if hit_count == record.alts().len() {
            failed_sites += 1;
            contained_sites += 1;
            if vrecord.filter == "PASS" {
                vrecord.filter = "contained".to_string();
            }
        }

        passed_sites += 1;

        if vrecord.filter == "PASS" {
            pass_bed
                .write(format!("{}\t{}\t{}\n", record.chr(), record.pos(), record.end()).as_bytes())
                .expect("failed to write bed");
        } else {
            fail_bed
                .write(format!("{}\t{}\t{}\n", record.chr(), record.pos(), record.end()).as_bytes())
                .expect("failed to write bed");
        }

        vcf.write(format!("{}\n", vrecord).as_bytes())
            .expect("Failure to write record to vcf");
    }

    info!(
        "passed:{} failed:{} contained:{} mismatch:{} ",
        passed_sites, failed_sites, contained_sites, gt_missmatch
    );
}

#[cfg(test)]
mod tests {
    use std::vec;

    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_parse_mash() {
        let input = "chr10:100771810:100771812:CT:C  1|0:1|0:0|0:1|0:1|0:1|0:0|0:0|0:1|0:0|0:        dragen-ilmn".to_string();
        let results = parse_record(&input);

        let mut expected = MashRecord {
            chrom: "chr10".to_string(),
            start: "100771810".to_string(),
            end: "100771812".to_string(),
            ref_allele: "CT".to_string(),
            alt_alleles: vec!["C".to_string()],
            sources_lookup: HashSet::new(),
            genotype_sets: HashSet::new(),
        };
        expected.sources_lookup.insert("dragen-ilmn".to_string());

        assert_eq!(results.ref_allele, expected.ref_allele);
        assert_eq!(results.alt_alleles, expected.alt_alleles);
        assert_eq!(results.chrom, expected.chrom);
    }
    #[test]
    fn test_gt_mismatch() {
        let input = "chr10:101285645:101285668:TGGGGAGCTCCTGGGCTCAGTGG:T     1|1:1|1:0|1:1|1:1|1:1|1:0|1:0|1:1|0:1|1:,1|1:1|1:1|1:1|1:1|1:1|1:1|1:1|1:1|1:1|1:       dragen-ilmn,dv-hifi".to_string();
        let results = parse_record(&input);
        assert_eq!(results.genotype_conflict(), true);
    }
    #[test]
    fn test_alt_count() {
        let input = "chr10:101390764:101390769:AAAAG:A,G     0|2:0|2:1|2:0|2:0|0:0|0:1|0:1|2:0|1:2|0:        dv-hifi".to_string();
        let results = parse_record(&input);
        assert_eq!(results.alts().len(), 2);
    }
}
