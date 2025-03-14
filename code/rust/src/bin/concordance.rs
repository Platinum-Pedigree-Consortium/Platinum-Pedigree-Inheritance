use std::collections::{HashMap, HashSet};
use std::io::Read as IoRead;
use std::mem::swap;
use std::path::Path;
use std::str;

use clap::Parser;

use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Format, Header, Read, Reader, Writer};

use concordance::iht;
use concordance::iht::InheritanceBlock;
/// Filter variants based on expected haplotypes
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// CSV containing the inheritance vectors (optionally gzipped)
    #[arg(short, long)]
    inheritance: String,

    /// VCF file
    #[arg(short, long)]
    vcf: String,

    #[arg(short, long)]
    prefix: String,

    /// Sample ID of the father
    #[arg(short, long)]
    father: String,

    /// Skip concordance filtering, but keep all other filters
    #[arg(short, long, default_value_t = false)]
    skip: bool,

    /// Sample ID of the mother
    #[arg(short, long)]
    mother: String,

    /// Minimum quality score
    #[arg(short, long, default_value_t = 20.0)]
    qual: f32,
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

fn header_to_idx(h: &HeaderView) -> HashMap<String, usize> {
    let mut lookup = HashMap::new();
    for (i, mut x) in (*h).samples().into_iter().enumerate() {
        let mut s = String::new();
        x.read_to_string(&mut s).expect("converting name to idx"); // Read sample name in to `s`
        println!("idx:{} sample:{}", i, s); // output sample name
        lookup.insert(s, i);
    }
    return lookup;
}

fn concordant(parents: [[u32; 4]; 4], block: &InheritanceBlock, genos: Vec<u32>) -> i8 {
    for (i, c) in parents.iter().enumerate() {
        let mut genovec: Vec<u32> = Vec::new();
        for p in &block.parental_hap {
            let mut first_allele = c[allele_conversion(p.chars().nth(0).unwrap())];
            let mut second_allele = c[allele_conversion(p.chars().nth(1).unwrap())];

            if first_allele > second_allele {
                swap(&mut first_allele, &mut second_allele);
            }
            genovec.push(first_allele);
            genovec.push(second_allele);
        }
        if genovec == genos {
            return i as i8;
        }
    }
    return -1;
}

fn build_phased_haplotypes(
    parents: &[u32; 4],
    block: &InheritanceBlock,
) -> HashMap<String, [u32; 2]> {
    let mut sample_to_alleles: HashMap<String, [u32; 2]> = HashMap::new();

    let sample_names = block.samples.clone();

    for (i, p) in block.parental_hap.iter().enumerate() {
        sample_to_alleles.insert(
            sample_names[i].clone(),
            [
                parents[allele_conversion((*p).chars().nth(0).unwrap())],
                parents[allele_conversion((*p).chars().nth(1).unwrap())],
            ],
        );
    }

    return sample_to_alleles;
}

fn reorder_alleles(alleles: Vec<&[u8]>) -> (HashMap<usize, usize>, Vec<&[u8]>) {
    let mut alts: Vec<(&str, usize)> = Vec::new();
    let mut allele_index_lookup = HashMap::new();
    let mut new_alleles: Vec<&[u8]> = Vec::new();
    new_alleles.push(alleles[0]);

    for (idx, alt) in alleles.iter().enumerate() {
        if idx == 0 {
            continue;
        }
        alts.push((str::from_utf8(alt).unwrap(), idx));
    }
    alts.sort_by_key(|x| x.0);

    // inserting reference (which is always zero)
    allele_index_lookup.insert(0, 0);

    for (idx, tup) in alts.iter().enumerate() {
        new_alleles.push(tup.0.as_bytes());
        allele_index_lookup.insert(tup.1, idx + 1);
    }
    return (allele_index_lookup, new_alleles);
}

fn main() {
    let args = Args::parse();
    let mut bcf = Reader::from_path(args.vcf).expect("Error opening vcf file.");
    let header = bcf.header().clone();
    let wheader: Header = Header::from_template(&header);
    let mut inheritance = iht::parse_inht(args.inheritance);
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
    let mut all_het: i64 = 0;
    let mut all_ref: i64 = 0;
    let mut no_con: i64 = 0;
    let mut nocall_geno: i64 = 0;
    let mut lowq: i64 = 0;
    let mut fail_counts: HashMap<String, i64> = HashMap::new();

    for (_i, record_result) in bcf.records().enumerate() {
        let record = record_result.expect("Fail to read record");
        let mut newrecord = record.clone();

        let alleles = record.alleles();

        let reordered_alleles = reorder_alleles(alleles);

        let gts = record.genotypes().expect("Error reading genotypes");
        let chr = std::str::from_utf8(header.rid2name(record.rid().unwrap()).expect("Invalid rid"))
            .unwrap();
        let mut current_block_idx: usize = 0;

        let mut block = iht::get_iht_block(
            &mut inheritance,
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

        let mut genovec: Vec<u32> = Vec::new();
        let mut failed_site = false;
        let mut alt_count: usize = 0;
        let mut het_count: usize = 0;

        let samples = block.as_mut().unwrap().samples.clone();

        let mut failed_vcf_record_processing = |current_block: &mut InheritanceBlock| {
            failed += 1;
            current_block.failing_count += 1;
            *fail_counts.entry("".to_string()).or_insert(0) += 1;
            outfailvcf.write(&record).unwrap();
        };

        let qual = record.qual();
        if qual < args.qual {
            failed_vcf_record_processing(block.as_mut().unwrap());
            lowq += 1;
            continue;
        }

        let mut no_missing = true;

        // samples are ordered by inheritance vector header
        for s in samples.iter() {
            let gt = gts.get(ped_idx_lookup[s]);

            // We are only working with diploid regions
            if gt.len() != 2 {
                failed_site = true;
                break;
            }

            let mut first_allele = gt[0];
            let mut second_allele = gt[1];

            if first_allele == GenotypeAllele::UnphasedMissing
                || first_allele == GenotypeAllele::PhasedMissing
                || second_allele == GenotypeAllele::UnphasedMissing
                || second_allele == GenotypeAllele::PhasedMissing
            {
                nocall_geno += 1;
                failed_site = true;
                no_missing = false;
                break;
            }
            // this will blow up unless the previous if statement is present
            if first_allele.index().unwrap() > 0 || second_allele.index().unwrap() > 0 {
                alt_count += 1;
            }

            if first_allele.index().unwrap() != second_allele.index().unwrap() {
                het_count += 1;
            }

            // ensuring that the alleles are ordered allowing us to test if we can get the right configuration.
            if first_allele.index() > second_allele.index() {
                swap(&mut first_allele, &mut second_allele);
            }
            genovec.push(first_allele.index().unwrap());
            genovec.push(second_allele.index().unwrap());
        }

        // if the alt count is zero, we have a hom-ref site.
        if alt_count == 0 && no_missing {
            all_ref += 1;
            failed_site = true;
        }

        if failed_site {
            failed_vcf_record_processing(block.as_mut().unwrap());
            continue;
        }

        let mother_gt = gts.get(ped_idx_lookup[&args.mother]);
        let father_gt = gts.get(ped_idx_lookup[&args.father]);

        let mut parent_allele_count: HashSet<u32> = HashSet::new();

        let configurations = [
            [
                father_gt[0].index().unwrap(),
                father_gt[1].index().unwrap(),
                mother_gt[0].index().unwrap(),
                mother_gt[1].index().unwrap(),
            ],
            [
                father_gt[1].index().unwrap(),
                father_gt[0].index().unwrap(),
                mother_gt[0].index().unwrap(),
                mother_gt[1].index().unwrap(),
            ],
            [
                father_gt[0].index().unwrap(),
                father_gt[1].index().unwrap(),
                mother_gt[1].index().unwrap(),
                mother_gt[0].index().unwrap(),
            ],
            [
                father_gt[1].index().unwrap(),
                father_gt[0].index().unwrap(),
                mother_gt[1].index().unwrap(),
                mother_gt[0].index().unwrap(),
            ],
        ];

        parent_allele_count.insert(father_gt[0].index().unwrap());
        parent_allele_count.insert(father_gt[1].index().unwrap());
        parent_allele_count.insert(mother_gt[0].index().unwrap());
        parent_allele_count.insert(mother_gt[1].index().unwrap());

        // if everything is het we check that the parent alleles are the same too. We do this because if there are more than two alleles in the parents we can distinguish inheritance.
        if het_count == samples.len() && parent_allele_count.len() == 2 {
            all_het += 1;
            failed_site = true;
        }

        let con = concordant(configurations, &block.as_ref().unwrap(), genovec);

        // unable to find a genotype configuration that matches the inheritance vector
        if con == -1 && !args.skip {
            failed_site = true;
            no_con += 1;
        }

        if failed_site {
            failed_vcf_record_processing(block.as_mut().unwrap());
            continue;
        }

        if !args.skip {
            let mut new_gts: Vec<GenotypeAllele> = Vec::new();
            let sample_count = usize::try_from(record.sample_count()).unwrap();

            let phased =
                build_phased_haplotypes(&configurations[con as usize], &block.as_ref().unwrap());

            for sample_index in 0..sample_count {
                let sample_name =
                    String::from_utf8(header.samples()[sample_index].to_vec()).unwrap();
                let gt = gts.get(sample_index);

                if phased.contains_key(&sample_name) {
                    let pg = phased.get(&sample_name).unwrap();

                    let mut a0: usize = pg[0].try_into().unwrap();
                    let mut a1: usize = pg[1].try_into().unwrap();

                    a0 = *reordered_alleles.0.get(&a0).unwrap();
                    a1 = *reordered_alleles.0.get(&a1).unwrap();

                    new_gts.push(GenotypeAllele::Unphased(a0 as i32));
                    new_gts.push(GenotypeAllele::Phased(a1 as i32));
                } else {
                    for g in gt.iter() {
                        new_gts.push(g.clone());
                    }
                }
            }

            newrecord.push_genotypes(&new_gts).unwrap();
            newrecord.set_alleles(&reordered_alleles.1).unwrap();

            passed += 1;
            outvcf.write(&newrecord).unwrap();
            inheritance[current_block_idx].passing_count += 1;
        } else {
            outvcf.write(&record).unwrap();
        }
    }

    for (ov, c) in &fail_counts {
        println!("failed_reason\t{}\t{}", ov, c);
    }

    for b in &inheritance {
        println!(
            "block_stats\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            args.prefix,
            b.chrom,
            b.start,
            b.end,
            b.passing_count,
            b.failing_count,
            b.passing_count + b.failing_count,
            b.inherited_haps.len(),
        );
    }

    let con_rate = (passed as f64) / (passed + no_con) as f64;

    println!(
        "not in block: {} passed: {} failed: {} all-het: {} all-ref: {} non-concordant: {} nocall-geno: {} lq: {} concordance-rate: {:.3} ",
        block_fail, passed, failed, all_het, all_ref, no_con, nocall_geno, lowq, con_rate
    );
}
