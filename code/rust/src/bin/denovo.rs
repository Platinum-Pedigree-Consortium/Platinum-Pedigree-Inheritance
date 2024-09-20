use clap::Parser;
use core::f32;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Header, Read, Reader};
use std::collections::HashMap;
use std::fs::read_to_string;

/// A tool to filter any pedigree for denovo variants
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// CSV containing pedigree information
    #[arg(short, long)]
    ped: String,

    /// VCF file
    #[arg(short, long)]
    vcf: String,

    /// Min GQ
    #[arg(short, long, default_value_t = 10)]
    gq: i32,

    /// Min DP
    #[arg(short, long, default_value_t = 5)]
    dp: i32,

    /// Max unknown transmission events. e.g. a genotyping is missing
    #[arg(short, long, default_value_t = 3)]
    unknown: i32,

    /// Minimum ref allele fraction in either parent
    #[arg(short, long, default_value_t = 0.9)]
    min_parental_ref_frac: f32,

    /// most ancestral individuals that have parents/children.
    #[arg(short, long)]
    founders: String,
}

#[derive(Debug, Clone)]
enum Sex {
    Male,
    Female,
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Geno {
    HomRef,
    Het,
    HomAlt,
    NoCall,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Iht {
    Unknown,
    CanDenovo,
    IhtDenovo,
    Violation,
    Valid,
}

impl Default for Iht {
    fn default() -> Self {
        Iht::Unknown
    }
}

#[derive(Debug, Clone)]
struct Member {
    family: String,
    id: String,
    father: String,
    mother: String,
    pub children: Vec<String>,
    sex: Sex,
    genotype: (GenotypeAllele, GenotypeAllele),
    geno: Geno,
    gq: i32,
    ref_ratio: f32,
}

fn parse_ped(filename: &str) -> HashMap<String, Member> {
    let mut result: HashMap<String, Member> = HashMap::new();
    let mut children: HashMap<String, Vec<String>> = HashMap::new();

    for line in read_to_string(filename).unwrap().lines() {
        let parts: Vec<String> = line.split("\t").map(str::to_string).collect();

        let mut sample_sex = Sex::Unknown;

        if parts[4] == "1".to_string() {
            sample_sex = Sex::Male;
        }

        result.insert(
            parts[1].to_string(),
            Member {
                family: parts[0].to_string(),
                id: parts[1].to_string(),
                father: parts[2].to_string(),
                mother: parts[3].to_string(),
                children: Vec::new(),
                sex: sample_sex,
                genotype: (
                    GenotypeAllele::UnphasedMissing,
                    GenotypeAllele::UnphasedMissing,
                ),
                geno: Geno::NoCall,
                gq: 0,
                ref_ratio: 1.0,
            },
        );

        if parts[2].to_string() == "NA".to_string() || parts[3].to_string() == "NA".to_string() {
            continue;
        }

        children
            .entry(parts[2].to_string())
            .or_default()
            .push(parts[1].to_string());
        children
            .entry(parts[3].to_string())
            .or_default()
            .push(parts[1].to_string());
    }

    for (k, v) in children {
        result.get_mut(&k).unwrap().children = v;
    }

    result
}

fn transmitted(parent: &Member, allele: &GenotypeAllele) -> i32 {
    let mut r = 0;
    if parent.genotype.0 == *allele {
        r += 1;
    }
    if parent.genotype.1 == *allele {
        r += 1;
    }
    return r;
}

fn check_iht(child: &Member, mother: &Member, father: &Member) -> (String, Iht) {
    if child.geno == Geno::NoCall || mother.geno == Geno::NoCall || father.geno == Geno::NoCall {
        return (child.id.clone(), Iht::Unknown);
    }

    if transmitted(mother, &child.genotype.0) > 0 && transmitted(father, &child.genotype.1) > 0 {
        return (child.id.clone(), Iht::Valid);
    }

    if transmitted(mother, &child.genotype.1) > 0 && transmitted(father, &child.genotype.0) > 0 {
        return (child.id.clone(), Iht::Valid);
    }

    if mother.geno == Geno::HomRef && father.geno == Geno::HomRef && child.geno == Geno::Het {
        return (child.id.clone(), Iht::CanDenovo);
    }
    /*
    if mother.geno == Geno::HOMALT && father.geno == Geno::HOMALT && child.geno == Geno::HET {
        return (child.id.clone(), Iht::CAN_DENOVO);
    }
    */

    return (child.id.clone(), Iht::Violation);
}

fn follow_family_transmission(
    family: &HashMap<String, Member>,
    focal: &String,
    transmissions: &mut HashMap<String, Iht>,
) {
    let person = family.get(focal).unwrap();
    if person.father == "NA".to_string() || person.mother == "NA".to_string() {
        return;
    }

    let father = family.get(&person.father).unwrap();
    let mother = family.get(&person.mother).unwrap();

    let mut status = check_iht(person, mother, father);

    let mut allele_counts: HashMap<GenotypeAllele, i32> = HashMap::new();
    allele_counts.insert(person.genotype.0, 0);
    allele_counts.insert(person.genotype.1, 0);
    for childk in &person.children {
        let child_member = family.get(childk).unwrap();
        allele_counts
            .entry(child_member.genotype.0)
            .and_modify(|counter| *counter += 1);
        allele_counts
            .entry(child_member.genotype.1)
            .and_modify(|counter| *counter += 1);
    }

    let mut acounter = 0;

    if status.1 == Iht::CanDenovo {
        for (_a, c) in allele_counts {
            if c > 0 {
                acounter += 1;
            }
        }
    }

    if acounter > 1 {
        status.1 = Iht::IhtDenovo;
    }

    transmissions.insert(status.0, status.1);

    for child in &person.children {
        follow_family_transmission(family, child, transmissions);
    }
}

fn avg_gq(family: &HashMap<String, Member>) -> f32 {
    let mut gq_sum = 0;
    for (_k, v) in family {
        gq_sum += v.gq;
    }

    if gq_sum == 0 {
        return 0f32;
    } else {
        return gq_sum as f32 / family.len() as f32;
    }
}

/*
fn lq_count(family: &HashMap<String, Member>, gq: i32) -> i32 {
    let mut lwcount = 0;
    for (_k, v) in family {
        if v.gq < gq {
            lwcount += 1;
        }
    }
    lwcount
}
*/

fn type_count(trans: &HashMap<String, Iht>) -> HashMap<String, i32> {
    let mut cm: HashMap<String, i32> = HashMap::new();

    cm.insert("CAN_DENOVO".to_string(), 0);
    cm.insert("IHT_DENOVO".to_string(), 0);
    cm.insert("UNKNOWN".to_string(), 0);
    cm.insert("VIOLATION".to_string(), 0);
    cm.insert("Valid".to_string(), 0);

    for v in trans {
        match v.1 {
            Iht::CanDenovo => *cm.entry("CAN_DENOVO".to_string()).or_insert(0) += 1,
            Iht::IhtDenovo => *cm.entry("IHT_DENOVO".to_string()).or_insert(0) += 1,
            Iht::Unknown => *cm.entry("UNKNOWN".to_string()).or_insert(0) += 1,
            Iht::Violation => *cm.entry("VIOLATION".to_string()).or_insert(0) += 1,
            Iht::Valid => *cm.entry("Valid".to_string()).or_insert(0) += 1,
        }
    }
    cm
}

#[derive(Debug, Clone)]
struct SiteInfo {
    denovo: bool,
    denovo_sample: String,
    classification: Iht,
    unknown: i32,
    total_count: i32,
    violations: i32,
    Valid: i32,
    avg_gq: f32,
    min_ref_fraction: f32,
}

fn get_first_denovo(trans: &HashMap<String, Iht>) -> String {
    let mut result = "NA".to_string();
    for (k, v) in trans {
        if *v == Iht::CanDenovo || *v == Iht::IhtDenovo {
            result = k.clone();
        }
    }
    result
}

fn process_site(family: &HashMap<String, Member>, gq: i32, founders: &Vec<String>) -> SiteInfo {
    let mut trans: HashMap<String, Iht> = HashMap::new();
    for k in founders {
        follow_family_transmission(family, k, &mut trans);
    }

    let counts = type_count(&trans);
    let dsample = get_first_denovo(&trans);
    let mut denovo = false;
    let mut classification = Iht::Unknown;
    let mut min_ref_fraction: f32 = 1.0;

    if dsample != "NA" {
        let person = family.get(&dsample).unwrap();
        denovo = true;
        classification = trans.get(&dsample).unwrap().clone();
        min_ref_fraction = f32::min(
            family.get(&person.father).unwrap().ref_ratio,
            family.get(&person.mother).unwrap().ref_ratio,
        );
    }

    if *counts.get("VIOLATION").unwrap() > 0 {
        denovo = false;
    }

    SiteInfo {
        denovo: denovo,
        denovo_sample: get_first_denovo(&trans),
        classification: classification,
        unknown: *counts.get("UNKNOWN").unwrap(),
        total_count: trans.len() as i32,
        violations: *counts.get("VIOLATION").unwrap(),
        Valid: *counts.get("Valid").unwrap(),
        avg_gq: avg_gq(family),
        min_ref_fraction: min_ref_fraction,
    }
}

fn main() {
    let args = Args::parse();
    let mut bcf = Reader::from_path(args.vcf).expect("Error opening vcf file.");
    let header = bcf.header().clone();

    let founders: Vec<String> = args.founders.split(",").map(str::to_string).collect();

    let kindred = parse_ped(&args.ped);

    let sample_count = usize::try_from(header.sample_count()).expect("failure to get samples");
    let mut sample_lookup: HashMap<String, usize> = HashMap::new();

    for sample_index in 0..sample_count {
        let sample_name = String::from_utf8(header.samples()[sample_index].to_vec()).unwrap();
        sample_lookup.insert(sample_name, sample_index);
    }

    println!("#chr\tpos\talleles\tallele_count\tdenovo_sample\tclassification\tValid_count\tviolation_count\tunknown_count\ttotal_transmissions\taverage_site_gq\tmin_ref_fraction");

    for (_i, record_result) in bcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let chr = std::str::from_utf8(header.rid2name(record.rid().unwrap()).expect("InValid rid"))
            .unwrap();

        // Loading up genotypes into the family structure.
        let mut local_family = kindred.clone();
        let gts = record.genotypes().expect("Error reading genotypes");
        let gqs = record.format(b"GQ").integer();
        let dps = record.format(b"DP").integer();
        let ad_values = record.format(b"AD").integer().unwrap();

        for (sample, idx) in &sample_lookup {
            let sgt = gts.get(*idx);

            if !local_family.contains_key(sample) {
                continue;
            }

            let mut gq = match gqs {
                Ok(ref gqs) => gqs[*idx][0],
                Err(_) => 0i32,
            };

            let dp: i32 = match dps {
                Ok(ref dps) => dps[*idx][0],
                Err(_) => 0i32,
            };

            let sf = &mut local_family.get_mut(sample).unwrap();

            if sgt.len() != 2 {
                (*sf).geno = Geno::NoCall;
                continue;
            }
            (*sf).genotype.0 = sgt[0];
            (*sf).genotype.1 = sgt[1];
            if sgt[0] == GenotypeAllele::UnphasedMissing
                || sgt[1] == GenotypeAllele::UnphasedMissing
                || dp < args.dp
                || gq < args.gq
            {
                gq = 0;
                (*sf).geno = Geno::NoCall;
            } else if sgt[0] != sgt[1] {
                (*sf).geno = Geno::Het;
            } else if sgt[0] == GenotypeAllele::Unphased(0) {
                (*sf).geno = Geno::HomRef;
            } else {
                (*sf).geno = Geno::HomAlt;
            }
            (*sf).gq = gq;

            if (*sf).geno != Geno::NoCall {
                let ref_count: f32 = ad_values[*idx][0] as f32;
                let mut total_count: f32 = 0.0;
                for v in ad_values[*idx] {
                    total_count += *v as f32;
                }
                if ref_count == 0.0 {
                    (*sf).ref_ratio = 0.0;
                } else {
                    (*sf).ref_ratio = ref_count / total_count;
                }
            }
        }
        let mut allele_string = String::new();
        for allele in record.alleles() {
            for c in allele {
                allele_string.push(char::from(*c))
            }
            allele_string.push(',');
        }

        let site_res = process_site(&local_family, args.gq, &founders);
        if site_res.denovo
            && site_res.unknown < args.unknown
            && args.min_parental_ref_frac < site_res.min_ref_fraction
        {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{:?}\t{}\t{}\t{}\t{}\t{}\t{}",
                chr,
                record.pos() + 1,
                allele_string,
                record.allele_count(),
                site_res.denovo_sample,
                site_res.classification,
                site_res.Valid,
                site_res.violations,
                site_res.unknown,
                site_res.total_count,
                site_res.avg_gq,
                site_res.min_ref_fraction,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    // Note this useful idiom: importing names from outer (for mod tests) scope.
    use super::*;

    #[test]
    fn test_trio_het0() {
        let mom = Member {
            family: "NA".to_string(),
            id: "MOM".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Female,
            genotype: (GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)),
            geno: Geno::HomRef,
            gq: 0,
            ref_ratio: 1.0,
        };

        let dad = Member {
            family: "NA".to_string(),
            id: "DAD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)),
            geno: Geno::HomAlt,
            gq: 0,
            ref_ratio: 1.0,
        };

        let child = Member {
            family: "NA".to_string(),
            id: "CHILD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        assert_eq!(check_iht(&child, &mom, &dad).1, Iht::Valid);
    }

    #[test]
    fn test_trio_het1() {
        let mom = Member {
            family: "NA".to_string(),
            id: "MOM".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Female,
            genotype: (GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        let dad = Member {
            family: "NA".to_string(),
            id: "DAD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        let child = Member {
            family: "NA".to_string(),
            id: "CHILD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        assert_eq!(check_iht(&child, &mom, &dad).1, Iht::Valid);
    }
    #[test]
    fn test_trio_het2() {
        let mom = Member {
            family: "NA".to_string(),
            id: "MOM".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Female,
            genotype: (GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)),
            geno: Geno::HomAlt,
            gq: 0,
            ref_ratio: 1.0,
        };

        let dad = Member {
            family: "NA".to_string(),
            id: "DAD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        let child = Member {
            family: "NA".to_string(),
            id: "CHILD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(0)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        assert_eq!(check_iht(&child, &mom, &dad).1, Iht::Valid);
    }
    #[test]
    fn test_trio_het3() {
        let mom = Member {
            family: "NA".to_string(),
            id: "MOM".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Female,
            genotype: (GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        let dad = Member {
            family: "NA".to_string(),
            id: "DAD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)),
            geno: Geno::HomRef,
            gq: 0,
            ref_ratio: 1.0,
        };

        let child = Member {
            family: "NA".to_string(),
            id: "CHILD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(0)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        assert_eq!(check_iht(&child, &mom, &dad).1, Iht::Valid);
    }

    #[test]
    fn test_trio_denovo1() {
        let mom = Member {
            family: "NA".to_string(),
            id: "MOM".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Female,
            genotype: (GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)),
            geno: Geno::HomAlt,
            gq: 0,
            ref_ratio: 1.0,
        };

        let dad = Member {
            family: "NA".to_string(),
            id: "DAD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)),
            geno: Geno::HomAlt,
            gq: 0,
            ref_ratio: 1.0,
        };

        let child = Member {
            family: "NA".to_string(),
            id: "CHILD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(0)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        assert_eq!(check_iht(&child, &mom, &dad).1, Iht::Violation);
    }

    #[test]
    fn test_trio_denovo2() {
        let mom = Member {
            family: "NA".to_string(),
            id: "MOM".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Female,
            genotype: (GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)),
            geno: Geno::HomRef,
            gq: 0,
            ref_ratio: 1.0,
        };

        let dad = Member {
            family: "NA".to_string(),
            id: "DAD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(0)),
            geno: Geno::HomRef,
            gq: 0,
            ref_ratio: 1.0,
        };

        let child = Member {
            family: "NA".to_string(),
            id: "CHILD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(0)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        assert_eq!(check_iht(&child, &mom, &dad).1, Iht::CanDenovo);
    }

    #[test]
    fn test_trio_twoalt() {
        let mom = Member {
            family: "NA".to_string(),
            id: "MOM".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Female,
            genotype: (GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(2)),
            geno: Geno::HomAlt,
            gq: 0,
            ref_ratio: 1.0,
        };

        let dad = Member {
            family: "NA".to_string(),
            id: "DAD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(3), GenotypeAllele::Unphased(4)),
            geno: Geno::HomAlt,
            gq: 0,
            ref_ratio: 1.0,
        };

        let child = Member {
            family: "NA".to_string(),
            id: "CHILD".to_string(),
            father: "NA".to_string(),
            mother: "NA".to_string(),
            children: Vec::new(),
            sex: Sex::Male,
            genotype: (GenotypeAllele::Unphased(3), GenotypeAllele::Unphased(2)),
            geno: Geno::Het,
            gq: 0,
            ref_ratio: 1.0,
        };

        assert_eq!(check_iht(&child, &mom, &dad).1, Iht::Valid);
    }
}
