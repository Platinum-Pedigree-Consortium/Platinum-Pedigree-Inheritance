use std::collections::HashMap;
use std::io::Read as IoRead;
use std::path::Path;

use clap::Parser;
use csv::ReaderBuilder;
use edit_distance::edit_distance;
use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Format, Header, Read, Reader, Writer};

/// Filter variants based on expected haplotypes
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// CSV containing the inheritance vectors
    #[arg(short, long)]
    inheritance: String,

    /// VCF file
    #[arg(short, long)]
    vcf: String,

    #[arg(short, long)]
    prefix: String,
}

const POSSIBLE_PATTERNS: [[i32; 4]; 15] = [
    //[ 0, 0, 0, 0], // dad homref mom homref
    [0, 0, 0, 1], // dad homref mom het (1)
    [0, 0, 1, 0], // dad homref mom het (2)
    [0, 0, 1, 1], // dad homref mom homalt
    [0, 1, 0, 0], // dad het mom homref (1)
    [1, 0, 0, 0], // dad het mom homref (2)
    [1, 0, 1, 0], // dad het mom het (1)
    [0, 1, 0, 1], // dad het mom het (2)
    [1, 0, 0, 1], // dad het mom het (3)
    [0, 1, 1, 0], // dad het mom het (4)
    [0, 1, 1, 1], // dad het mom homalt (1)
    [1, 0, 1, 1], // dad het mom homalt (2)
    [1, 1, 0, 0], // dad is homalt mom is homref
    [1, 1, 0, 1], // dad is homalt mom is het (1)
    [1, 1, 1, 0], // dad is homalt mom is het (2)
    [1, 1, 1, 1], // dad is homalt mom is homalt
];

fn geno_conversion(geno: String) -> i32 {
    match geno.as_str() {
        "0/0" => 0,
        "0|0" => 0,
        "0/1" => 1,
        "1/0" => 1,
        "0|1" => 1,
        "1|0" => 1,
        "1/1" => 2,
        "1|1" => 2,
        _ => 3,
    }
}

fn allele_conversion(allele: char) -> usize {
    match allele {
        'A' => 0,
        'B' => 1,
        'C' => 2,
        'D' => 3,
        _ => 4,
    }
}

struct InheritanceBlock {
    chrom: String,
    start: i32,
    end: i32,
    passing_count: i32,
    failing_count: i32,
    samples: Vec<String>,
    sample_lookups: HashMap<String, usize>,
    parental_hap: Vec<String>,
    patterns: HashMap<String, Vec<[i32; 2]>>,
}
// returns a HashMap were the key is a simplified genotype string (0=homref,1=het,2=homalt), and the corresponding phase
fn okay_genotypes(iht: &InheritanceBlock) -> HashMap<String, Vec<[i32; 2]>> {
    let mut lookup: HashMap<String, Vec<[i32; 2]>> = HashMap::new();
    // possible genotypes
    for i in POSSIBLE_PATTERNS {
        let mut onevec: String = "".to_string();
        let mut phased_genotypes: Vec<[i32; 2]> = Vec::new();

        // possible haplotype configurations
        for j in &iht.parental_hap {
            let possible_geno = format!(
                "{}|{}",
                i[allele_conversion(j.chars().next().unwrap())],
                i[allele_conversion(j.chars().nth(1).unwrap())]
            );
            phased_genotypes.push([
                i[allele_conversion(j.chars().next().unwrap())],
                i[allele_conversion(j.chars().nth(1).unwrap())],
            ]);

            onevec = format!("{}{}", onevec, geno_conversion(possible_geno))
        }

        if lookup.contains_key(&onevec) {
            let seen_phase = lookup.get(&onevec).unwrap();
            for (g1, g2) in seen_phase.iter().zip(phased_genotypes.iter()) {
                if g1 != g2 {
                    println!(
                        "WARNING disagreement in phase: {} {:?} {:?}",
                        onevec, seen_phase, phased_genotypes
                    );
                }
            }
        } else {
            lookup.insert(onevec.clone(), phased_genotypes);
        }
    }
    lookup
}

impl std::fmt::Display for InheritanceBlock {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "inheritance-block: {} {} {} {:?} {:?} {:?}",
            self.chrom, self.start, self.end, self.samples, self.parental_hap, self.patterns
        )
    }
}

fn parse_inht(inht_fn: String) -> Vec<InheritanceBlock> {
    let mut inht_info = Vec::new();
    let mut reader = ReaderBuilder::new().from_path(inht_fn);

    let header = reader
        .as_mut()
        .expect("Error reading inheritance CSV header")
        .headers()
        .expect("Error reading inheritance CSV header")
        .clone();

    for record in reader.expect("Error reading inheritance CSV.").records() {
        let mut ihtblock = InheritanceBlock {
            chrom: record.as_ref().unwrap()[0].to_string().clone(),
            start: record.as_ref().unwrap()[1].parse::<i32>().unwrap(),
            end: record.as_ref().unwrap()[2].parse::<i32>().unwrap(),
            passing_count: 0,
            failing_count: 0,
            samples: Vec::new(),
            sample_lookups: HashMap::new(),
            parental_hap: Vec::new(),
            patterns: HashMap::new(),
        };

        let mut one_parent = false;

        for (sidx, i) in (3..header.len()).enumerate() {
            // Some inheritance vectors do not contain both parental marker, for example you might see `A` rather than `AB`.
            // We want to skip these sites as the expected genotypes are unknown.
            if record.as_ref().unwrap()[i].to_string().len() == 1 {
                one_parent = true;
            }
            // putting the sample names in the header into the inheritance block for easy lookup
            ihtblock.sample_lookups.insert(header[i].to_string(), sidx);
            ihtblock.samples.push(header[i].to_string());
            // geno
            ihtblock
                .parental_hap
                .push(record.as_ref().unwrap()[i].to_string());
        }
        if one_parent {
            println!("Warning skipping block missing both parents {}", ihtblock);
            continue;
        }
        ihtblock.patterns = okay_genotypes(&ihtblock);
        println!("{}", ihtblock);
        inht_info.push(ihtblock);
    }
    inht_info
}

fn header_to_idx(h: &HeaderView) -> HashMap<String, usize> {
    let mut lookup = HashMap::new();
    for (i, mut x) in (*h).samples().into_iter().enumerate() {
        let mut s = String::new();
        x.read_to_string(&mut s).expect("converting name to idx"); // Read sample name in to `s`
        println!("idx:{} sample:{}", i, s); // output sample name
        lookup.insert(s, i);
    }
    lookup
}

fn get_iht_block<'a>(
    ihts: &'a [InheritanceBlock],
    chr: &str,
    pos: i32,
    current: &mut usize,
) -> Option<&'a InheritanceBlock> {
    for item in ihts.iter().skip(*current) {
        if (item.chrom == chr) && (pos >= item.start) && (pos <= item.end) {
            return Some(item);
        }
        *current += 1;
    }
    None
}

fn main() {
    let args = Args::parse();
    let mut bcf = Reader::from_path(args.vcf).expect("Error opening vcf file.");
    let header = bcf.header().clone();
    let wheader: Header = Header::from_template(&header);
    let mut inheritance = parse_inht(args.inheritance);
    let ped_idx_lookup = header_to_idx(&header);
    let mut outvcf = Writer::from_path(
        Path::new(&format!("{}.vcf", args.prefix)),
        &wheader,
        true,
        Format::Vcf,
    )
    .unwrap();

    let mut outfailvcf = Writer::from_path(
        Path::new(&format!("{}.failedsites.vcf", args.prefix)),
        &wheader,
        true,
        Format::Vcf,
    )
    .unwrap();

    let mut block_fail: i64 = 0;
    let mut failed: i64 = 0;
    let mut passed: i64 = 0;
    let mut fail_counts: HashMap<String, i64> = HashMap::new();

    for record_result in bcf.records() {
        let mut record = record_result.expect("Fail to read record");
        // only looking at bi-allelic sites
        if record.alleles().len() != 2 {
            println!("Warning: skipping multi-allelic variant!");
            continue;
        }
        let gts = record.genotypes().expect("Error reading genotypes");
        let chr = std::str::from_utf8(header.rid2name(record.rid().unwrap()).expect("Invalid rid"))
            .unwrap();
        let mut current_block_idx: usize = 0;

        let block = get_iht_block(
            &inheritance,
            chr,
            record.pos().try_into().unwrap(),
            &mut current_block_idx,
        );
        match block {
            None => {
                println!(
                    "Warning: skipping variant as it is not in a block {} {}",
                    chr,
                    record.pos()
                );
                block_fail += 1;
                continue;
            }
            Some(_) => {}
        }
        let mut onevec: String = "".to_string();
        let mut genostring: String = "".to_string();

        for s in block.unwrap().samples.iter() {
            onevec = format!(
                "{}{}",
                onevec,
                geno_conversion(gts.get(ped_idx_lookup[s]).to_string())
            );
            genostring = format!("{},{}", genostring, gts.get(ped_idx_lookup[s]));
        }

        // this is a hom ref site for the individuals in the group of interest.
        let novar = String::from_utf8(vec![b'0'; onevec.len()]).unwrap();
        // this is a het site across the individuals in the group of interest.
        let all_het = String::from_utf8(vec![b'1'; onevec.len()]).unwrap();

        // no variants at this site or all hets which
        if (onevec == novar) || (onevec == all_het) {
            continue;
        }

        if !block.unwrap().patterns.contains_key(&onevec) {
            failed += 1;
            inheritance[current_block_idx].failing_count += 1;
            *fail_counts.entry(onevec.clone()).or_insert(0) += 1;

            outfailvcf.write(&record).unwrap();

            let mut oneoffcount = 0;
            let mut closestmatch: String = "".to_string();

            for (k, _v) in inheritance[current_block_idx].patterns.iter() {
                let ed = edit_distance(k, &onevec);
                if ed == 1 {
                    oneoffcount += 1;
                    closestmatch = k.to_string();
                }
            }
            if oneoffcount == 1 {
                let mut violator = "NA".to_string();

                for s in 0..onevec.len() {
                    if onevec.as_bytes()[s] != closestmatch.as_bytes()[s] {
                        violator = inheritance[current_block_idx].samples[s].clone();
                    }
                }

                let has_nocall = onevec.find('3');
                let nocall_violation = match has_nocall {
                    Some(_index) => true,
                    None => false,
                };

                let mut indel = false;
                let alleles = record.alleles();
                if alleles[0].len() != alleles[1].len() {
                    indel = true;
                }

                println!(
                    "chr:{} pos:{} indel:{} vec:{} closest:{} violator:{} nocall:{}",
                    chr,
                    record.pos(),
                    indel,
                    onevec,
                    closestmatch,
                    violator,
                    nocall_violation,
                );
            }

            continue;
        }

        let mut new_gts: Vec<GenotypeAllele> = Vec::new();
        let sample_count = usize::try_from(record.sample_count()).unwrap();
        let phased_genos = block.unwrap().patterns.get(&onevec).unwrap();

        for sample_index in 0..sample_count {
            let sample_name = String::from_utf8(header.samples()[sample_index].to_vec()).unwrap();
            let gt = gts.get(sample_index);

            if block.unwrap().sample_lookups.contains_key(&sample_name) {
                let gt_idx: usize = *block.unwrap().sample_lookups.get(&sample_name).unwrap();

                new_gts.push(GenotypeAllele::Phased(phased_genos[gt_idx][0]));
                new_gts.push(GenotypeAllele::Phased(phased_genos[gt_idx][1]));

                let new_geno: String =
                    format!("{}/{}", phased_genos[gt_idx][0], phased_genos[gt_idx][1]);

                /*
                println!(
                    "Setting genotypes: pos:{}, expected:{}, observed:{}, sample:{}, sample_index:{}, gt_idx:{}, samples:{:?}, genos:{:?}, geno_string{:?}",
                    record.pos(),
                    gt.to_string(),
                    new_geno,
                    sample_name,
                    sample_index,
                    gt_idx,
                    block.unwrap().samples,
                    phased_genos,
                    genostring,
                );
                */

                if geno_conversion(new_geno.clone()) != geno_conversion(gt.to_string()) {
                    panic!("FATAL: trying to set the wrong genotype");
                }
            } else {
                for g in gt.iter().cloned() {
                    new_gts.push(g);
                }
            }
        }

        record.push_genotypes(&new_gts).unwrap();

        passed += 1;
        outvcf.write(&record).unwrap();
        inheritance[current_block_idx].passing_count += 1;
    }

    for (ov, c) in &fail_counts {
        println!("failed_reason\t{}\t{}", ov, c);
    }

    for b in &inheritance {
        println!(
            "block_stats\t{}\t{}\t{}\t{}\t{}",
            b.chrom, b.start, b.end, b.passing_count, b.failing_count
        );
    }

    println!(
        "not in block: {} passed: {} failed: {} ",
        block_fail, passed, failed
    );
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_geno_hom() {
        assert_eq!(geno_conversion("0/0".to_string()), 0);
    }
    #[test]
    fn test_geno_hom_phase() {
        assert_eq!(geno_conversion("0|0".to_string()), 0);
    }
    #[test]
    fn test_geno_het1() {
        assert_eq!(geno_conversion("0/1".to_string()), 1);
    }
    #[test]
    fn test_geno_het2() {
        assert_eq!(geno_conversion("1/0".to_string()), 1);
    }
    #[test]
    fn test_geno_het3() {
        assert_eq!(geno_conversion("1|0".to_string()), 1);
    }
    #[test]
    fn test_geno_het4() {
        assert_eq!(geno_conversion("0|1".to_string()), 1);
    }
    #[test]
    fn test_geno_nc1() {
        assert_eq!(geno_conversion("./1".to_string()), 3);
    }
    #[test]
    fn test_geno_nc2() {
        assert_eq!(geno_conversion("./.".to_string()), 3);
    }
    #[test]
    fn test_geno_nc3() {
        assert_eq!(geno_conversion("./0".to_string()), 3);
    }
}
