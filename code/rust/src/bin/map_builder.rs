use clap::Parser;
use concordance::bed::BedRecord;
use concordance::iht::ChromType;
use concordance::iht::Iht;
use concordance::iht::IhtVec;
use concordance::utils::*;

use itertools::Itertools;
use log::{debug, error, info, warn, LevelFilter};

use std::collections::{HashMap, HashSet};
use std::fmt;

use std::fs::OpenOptions;
use std::io;

use std::io::Read as IoRead;
use std::io::Write;

use std::path::Path;
use std::process;
use std::str;

use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::Record;
use rust_htslib::bcf::{IndexedReader, Read};

use concordance::ped::{Family, Individual};
/// Build a pedigree haplotype map (inheritance vectors).
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = " A tool to map haplotypes through generations"
)]
struct Args {
    /// A PED file documenting familial relationships
    #[arg(short, long)]
    ped: String,

    /// VCF file
    #[arg(long)]
    vcf: String,

    /// Output prefix
    #[arg(long)]
    prefix: String,

    /// Minimum variant quality score
    #[arg(short, long, default_value_t = 20.0)]
    qual: f32,

    /// Minimum depth for all family members
    #[arg(short, long, default_value_t = 10)]
    depth: i32,

    /// Shortest allowed haplotype (by marker count)
    #[arg(short, long, default_value_t = 10)]
    run: usize,

    /// Verbosity
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
}

fn check_individuals_in_vcf(ped_path: &str, vcf_path: &str) -> io::Result<()> {
    let family = Family::parse_ped_file(ped_path)?;

    let vcf_samples = get_samples(vcf_path).unwrap();

    for id in family.get_individuals_ids() {
        if !vcf_samples.contains(&id) {
            error!(
                "Individual ID {} from PED file is not found in the VCF header.",
                id
            );
            process::exit(1);
        }
    }
    debug!("Cross-validated vcf and ped files");
    Ok(())
}

fn is_vcf_indexed(vcf_path: &str) -> io::Result<bool> {
    let index_path = format!("{}.tbi", vcf_path);
    if Path::new(&index_path).exists() {
        Ok(true)
    } else {
        warn!("VCF file {} is not indexed.", vcf_path);
        Ok(false)
    }
}

/// Extract chromosome names (contig names) from the `HeaderView` of a VCF file.
fn extract_chromosome_names(
    vcf_path: &str,
) -> Result<HashMap<String, u32>, Box<dyn std::error::Error>> {
    let bcf = IndexedReader::from_path(vcf_path).expect("Error opening vcf file.");
    let header = bcf.header();

    // Iterate through all header records and extract contig entries
    let contig_count = header.contig_count();
    let mut chromosomes = HashMap::new();

    for i in 0..contig_count {
        let name = header
            .rid2name(i)
            .map_err(|_| format!("Failed to retrieve contig name for index {}", i))?;
        let name_str = str::from_utf8(name)
            .map_err(|_| format!("Failed to parse contig name for index {}", i))?;
        chromosomes.insert(name_str.to_string(), i);
    }

    Ok(chromosomes)
}

fn get_samples(vcf_path: &str) -> io::Result<Vec<String>> {
    let bcf = IndexedReader::from_path(vcf_path).expect("Error opening vcf file.");
    let header = bcf.header();

    let samples: Vec<String> = header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect();
    debug!("Loaded samples from VCF header OK.");
    Ok(samples)
}

/// Reads all records in a specific chromosome from a VCF file and extracts genotype alleles.
/// Returns a vector where each element is a HashMap for a VCF record.
/// The key is the sample ID, and the value is a vector of `GenotypeAllele`.
fn parse_vcf(
    reader: &mut IndexedReader,
    chromosome: u32,
    chomname: &String,
    minqual: f32,
) -> io::Result<Vec<(BedRecord, HashMap<String, (i32, Vec<GenotypeAllele>)>)>> {
    let rv = reader.fetch(chromosome, 0, None);
    // Initialize a vector to store genotype data for each site
    let mut records_genotype_map: Vec<(BedRecord, HashMap<String, (i32, Vec<GenotypeAllele>)>)> =
        Vec::new();

    match rv {
        Err(_) => return Ok(records_genotype_map),
        _ => {}
    }

    // Extract sample names from the VCF header
    let header = reader.header();
    let samples: Vec<String> = header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect();

    // Iterate through all records in the specified chromosome
    for result in reader.records() {
        let record = result.unwrap();

        if record.qual() < minqual {
            continue;
        }

        if is_indel(&record) {
            continue;
        }

        let depths = get_sample_depths(&record, &samples).unwrap();

        let genotypes = record.genotypes().expect("Failed to retrieve genotypes");

        let mut record_map: HashMap<String, (i32, Vec<GenotypeAllele>)> = HashMap::new();

        // Iterate through samples and extract genotype alleles
        for (i, sample) in samples.iter().enumerate() {
            let geno: rust_htslib::bcf::record::Genotype = genotypes.get(i);
            let alleles: Vec<GenotypeAllele> = geno.iter().cloned().collect(); // Clone each allele
            record_map.insert(sample.clone(), (*depths.get(sample).unwrap(), alleles));
        }
        records_genotype_map.push((
            BedRecord {
                chrom: chomname.clone(),
                start: record.pos(),
                end: record.pos(),
            },
            record_map,
        ));
    }

    Ok(records_genotype_map)
}

/// removes samples not in the keys list
fn remove_unused_samples(
    map: &mut HashMap<String, (i32, Vec<GenotypeAllele>)>,
    keys: &Vec<String>,
) {
    let key_set: std::collections::HashSet<_> = keys.iter().cloned().collect();
    map.retain(|key, _| key_set.contains(key));
}

/// Checks if any sample has missing alleles or a sequencing depth of zero.
///
/// # Arguments
/// - `genotype_data`: A `HashMap` mapping sample IDs to their depth and genotype alleles.
/// - `min depth` : A i32
/// # Returns
/// - `true` if any sample has a depth of less than 10 or missing.
/// - `false` if all samples have valid depths and alleles.
fn has_missing_alleles(
    genotype_data: &HashMap<String, (i32, Vec<GenotypeAllele>)>,
    min: i32,
) -> bool {
    genotype_data.values().any(|(depth, allele_vec)| {
        *depth < min || allele_vec.iter().any(|allele| allele.index().is_none())
    })
}

/// Identifies unique alleles for each founder and finds which children inherit them.
///
/// # Arguments
/// - `genotype_data`: A `HashMap` mapping individual IDs to `(i32, Vec<GenotypeAllele>)`.
/// - `family`: A reference to the `Family` struct.
///
/// # Returns
/// - A tuple where:
///   - `bool`: Indicates if a unique parent marker exists.
///   - `HashMap<String, (HashSet<GenotypeAllele>, Vec<String>)>`:
///     - Key: Founder ID.
///     - Value: `(Unique alleles, List of children with those alleles)`.
fn find_founder_unique_alleles(
    genotype_data: &HashMap<String, (i32, Vec<GenotypeAllele>)>,
    family: &Family,
) -> (
    bool,
    HashMap<String, (HashSet<GenotypeAllele>, Vec<String>)>,
) {
    let founders = family.founders(); // Retrieve founders from Family
    let children = family.offspring(); // Retrieve children from Family

    let mut founder_alleles: HashMap<String, HashSet<GenotypeAllele>> = HashMap::new();

    // Step 1: Collect alleles for each founder
    for founder in &founders {
        if let Some((_, alleles)) = genotype_data.get(&founder.id()) {
            founder_alleles
                .entry(founder.id())
                .or_insert_with(HashSet::new)
                .extend(alleles.iter().cloned());
        }
    }

    // Step 2: Identify unique alleles for each founder
    let mut unique_alleles_map: HashMap<String, HashSet<GenotypeAllele>> = founder_alleles.clone();

    for (founder, alleles) in &mut unique_alleles_map {
        // Remove alleles that are present in other founders
        for (other_founder, other_alleles) in &founder_alleles {
            if founder != other_founder {
                alleles.retain(|allele| !other_alleles.contains(allele));
            }
        }
    }

    // Step 3: Filter out founders without unique alleles
    unique_alleles_map.retain(|_, alleles| !alleles.is_empty());

    // If there's not exactly one unique founder, return empty results
    if unique_alleles_map.len() != 1 {
        return (false, HashMap::new());
    }

    let unique_founder_id = unique_alleles_map.keys().next().unwrap().clone();
    let unique_alleles = unique_alleles_map.get(&unique_founder_id).unwrap();

    let mut children_with_unique_alleles = Vec::new();

    // Step 4: Identify children inheriting unique alleles
    for child in &children {
        if let Some((_, child_alleles)) = genotype_data.get(&child.id()) {
            let has_unique_allele = child_alleles.iter().any(|a| unique_alleles.contains(a));

            if has_unique_allele {
                let father_id = child.get_father_id().unwrap_or_default();
                let mother_id = child.get_mother_id().unwrap_or_default();

                // Check if the unique allele is actually present in one of the parents
                let father_has_allele = genotype_data.get(&father_id).map_or(false, |(_, fa)| {
                    fa.iter().any(|a| unique_alleles.contains(a))
                });

                let mother_has_allele = genotype_data.get(&mother_id).map_or(false, |(_, ma)| {
                    ma.iter().any(|a| unique_alleles.contains(a))
                });

                // If neither parent has the unique allele, return empty
                if !father_has_allele && !mother_has_allele {
                    return (false, HashMap::new());
                }

                children_with_unique_alleles.push(child.id());
            }
        }
    }

    (
        !children_with_unique_alleles.is_empty(),
        HashMap::from([(
            unique_founder_id,
            (unique_alleles.clone(), children_with_unique_alleles),
        )]),
    )
}

fn load_up(
    fam: &Family,
    iht: &mut Iht,
    markers: &HashMap<String, (HashSet<GenotypeAllele>, Vec<String>)>,
    zygosity: &ChromType,
) {
    for f in markers {
        let founder_alleles = iht.get_alleles(f.0).unwrap();
        let founder_info = fam.get_individual(f.0).unwrap();

        // By default, use the first allele; if on ChrX and the founder is male, flip to the second allele.
        let mut founder_allele = founder_alleles.0;
        if *zygosity == ChromType::ChrX && founder_info.get_sex().unwrap() == 1 {
            founder_allele = founder_alleles.1;
        }

        for o in &f.1 .1 {
            let child_alleles = iht.children.get(o).unwrap();

            // This logic deals with allele inheritance across generations.
            let child = fam.get_individual(o).unwrap();
            let father = fam.get_individual(&child.get_father_id().unwrap()).unwrap();

            if f.1 .1.contains(&father.id()) || *f.0 == father.id() {
                // paternal allele goes to the left
                iht.children
                    .insert((*o).clone(), (founder_allele, child_alleles.1));
            } else {
                // maternal allele goes to the right
                iht.children
                    .insert((*o).clone(), (child_alleles.1, founder_allele));
            }
        }
    }

    // After processing markers, update children genotypes for ChrX:
    // If the ChromType is ChrX and the child is male, set the first genotype entry to "."
    if *zygosity == ChromType::ChrX {
        for (child_id, genotype) in iht.children.iter_mut() {
            let child = fam.get_individual(child_id).unwrap();
            if child.get_sex().unwrap() == 1 {
                genotype.0 = '.';
            }
        }
    }
}

/// Applies the threshold filtering over a mutable vector of IhtVecs.
fn apply_threshold_to_vec(iht_vecs: &mut Vec<IhtVec>, threshold: usize) {
    let mut prev: Option<&IhtVec> = None;
    let mut iter = iht_vecs.iter_mut().peekable();

    while let Some(curr) = iter.next() {
        let next = iter.peek().map(|x| &**x);
        apply_threshold(curr, threshold, prev, next);
        prev = Some(curr);
    }
}

/// Replaces characters with '?' if their count is below the threshold,
/// unless they match both the previous and next IhtVec.
fn apply_threshold(
    iht_vec: &mut IhtVec,
    threshold: usize,
    prev: Option<&IhtVec>,
    next: Option<&IhtVec>,
) {
    for (child, char_counts) in &iht_vec.non_missing_counts {
        if let Some((hap_a, hap_b)) = iht_vec.iht.children.get_mut(child) {
            let prev_match = prev.and_then(|p| p.iht.children.get(child)).cloned();
            let next_match = next.and_then(|n| n.iht.children.get(child)).cloned();

            let should_keep = |hap: char| {
                prev_match.map_or(false, |(pa, pb)| pa == hap || pb == hap)
                    && next_match.map_or(false, |(na, nb)| na == hap || nb == hap)
            };

            if char_counts.get(hap_a).copied().unwrap_or(0) < threshold && !should_keep(*hap_a) {
                *hap_a = '?';
            }
            if char_counts.get(hap_b).copied().unwrap_or(0) < threshold && !should_keep(*hap_b) {
                *hap_b = '?';
            }
        }
    }
}

fn collapse_identical_iht(data: Vec<IhtVec>) -> Vec<IhtVec> {
    let mut collapsed = Vec::new();
    let mut iter = data.into_iter().peekable();

    while let Some(current) = iter.next() {
        let start = current.bed.start;
        let mut end = current.bed.end;
        let mut count = current.count;
        let mut merged_founders = current.iht.founders.clone();
        let mut merged_children = current.iht.children.clone();
        let mut non_missing_counts = current.non_missing_counts;

        while let Some(next) = iter.peek() {
            if next.bed.chrom == current.bed.chrom
                && can_merge_families(
                    &merged_founders,
                    &next.iht.founders,
                    &merged_children,
                    &next.iht.children,
                )
            {
                end = next.bed.end;
                count += next.count;
                merge_family_maps(&mut merged_founders, &next.iht.founders);
                merge_family_maps(&mut merged_children, &next.iht.children);
                merge_non_missing_counts(&mut non_missing_counts, &next.non_missing_counts);
                iter.next();
            } else {
                break;
            }
        }

        collapsed.push(IhtVec {
            bed: BedRecord {
                chrom: current.bed.chrom.clone(),
                start,
                end,
            },
            count,
            iht: Iht {
                founders: merged_founders,
                children: merged_children,
            },
            non_missing_counts,
        });
    }

    collapsed
}

fn count_non_missing(
    founders: &HashMap<String, (char, char)>,
    children: &HashMap<String, (char, char)>,
) -> HashMap<String, HashMap<char, usize>> {
    let mut counts: HashMap<String, HashMap<char, usize>> = HashMap::new();

    for (key, &(c1, c2)) in founders.iter().chain(children.iter()) {
        let entry = counts.entry(key.clone()).or_insert_with(HashMap::new);
        if c1 != '?' {
            *entry.entry(c1).or_insert(0) += 1;
        }
        if c2 != '?' {
            *entry.entry(c2).or_insert(0) += 1;
        }
    }
    counts
}

fn merge_non_missing_counts(
    map1: &mut HashMap<String, HashMap<char, usize>>,
    map2: &HashMap<String, HashMap<char, usize>>,
) {
    for (key, char_counts) in map2.iter() {
        let entry = map1.entry(key.clone()).or_insert_with(HashMap::new);
        for (char_key, count) in char_counts.iter() {
            *entry.entry(*char_key).or_insert(0) += count;
        }
    }
}

fn can_merge_families(
    founders1: &HashMap<String, (char, char)>,
    founders2: &HashMap<String, (char, char)>,
    children1: &HashMap<String, (char, char)>,
    children2: &HashMap<String, (char, char)>,
) -> bool {
    founders1.iter().all(|(key, &(c1, c2))| {
        founders2.get(key).map_or(true, |&(d1, d2)| {
            (c1 == '?' || d1 == '?' || c1 == d1) && (c2 == '?' || d2 == '?' || c2 == d2)
        })
    }) && children1.iter().all(|(key, &(c1, c2))| {
        children2.get(key).map_or(true, |&(d1, d2)| {
            (c1 == '?' || d1 == '?' || c1 == d1) && (c2 == '?' || d2 == '?' || c2 == d2)
        })
    })
}

fn merge_family_maps(
    map1: &mut HashMap<String, (char, char)>,
    map2: &HashMap<String, (char, char)>,
) {
    for (key, &(c1, c2)) in map2.iter() {
        map1.entry(key.clone())
            .and_modify(|(m1, m2)| {
                if *m1 == '?' && c1 != '?' {
                    *m1 = c1;
                }
                if *m2 == '?' && c2 != '?' {
                    *m2 = c2;
                }
            })
            .or_insert((c1, c2));
    }
}

fn is_indel(record: &Record) -> bool {
    // Get reference allele
    let alleles = record.alleles();

    if alleles.len() < 2 {
        return false; // No alternate allele
    }

    let ref_len = alleles[0].len();
    for alt in &alleles[1..] {
        // Check all alternative alleles
        if alt.len() != ref_len {
            return true; // Length mismatch = indel
        }
    }

    false
}

/// Converts a HashMap<String, (HashSet<GenotypeAllele>, Vec<String>)> to a sorted string representation.
fn marker_to_string(map: &HashMap<String, (HashSet<GenotypeAllele>, Vec<String>)>) -> String {
    let mut entries: Vec<String> = map
        .iter()
        .map(|(key, (alleles, vec))| {
            let mut sorted_vec = vec.clone();
            sorted_vec.sort();
            format!(
                "{} {:?} {}",
                key,
                alleles.iter().collect::<Vec<_>>().get(0).unwrap(),
                sorted_vec.join(",")
            )
        })
        .collect();

    entries.sort();
    entries.join(", ")
}

/// Fills in '?' values when matching entries exist before and after, iterating until no more changes occur.
fn fill_missing_values(iht_vecs: &mut Vec<IhtVec>) {
    if iht_vecs.len() < 3 {
        return;
    }

    let mut changes_made;
    loop {
        changes_made = false;
        let mut to_update = Vec::new();

        for i in 0..iht_vecs.len() {
            for (child, (hap_a, hap_b)) in &iht_vecs[i].iht.children {
                if *hap_a == '?' || *hap_b == '?' {
                    let mut matching_hap_a = None;
                    let mut matching_hap_b = None;

                    for j in (0..i).rev() {
                        if let Some((prev_a, prev_b)) = iht_vecs[j].iht.children.get(child) {
                            if matching_hap_a.is_none() && *prev_a != '?' {
                                matching_hap_a = Some(*prev_a);
                            }
                            if matching_hap_b.is_none() && *prev_b != '?' {
                                matching_hap_b = Some(*prev_b);
                            }
                            if matching_hap_a.is_some() && matching_hap_b.is_some() {
                                break;
                            }
                        }
                    }

                    for j in i + 1..iht_vecs.len() {
                        if let Some((next_a, next_b)) = iht_vecs[j].iht.children.get(child) {
                            if matching_hap_a.is_none() && *next_a != '?' {
                                matching_hap_a = Some(*next_a);
                            }
                            if matching_hap_b.is_none() && *next_b != '?' {
                                matching_hap_b = Some(*next_b);
                            }
                            if matching_hap_a.is_some() && matching_hap_b.is_some() {
                                break;
                            }
                        }
                    }

                    let mut new_hap_a = *hap_a;
                    let mut new_hap_b = *hap_b;

                    if *hap_a == '?' && matching_hap_a.is_some() {
                        new_hap_a = matching_hap_a.unwrap();
                        changes_made = true;
                    }
                    if *hap_b == '?' && matching_hap_b.is_some() {
                        new_hap_b = matching_hap_b.unwrap();
                        changes_made = true;
                    }

                    if changes_made {
                        to_update.push((i, child.clone(), new_hap_a, new_hap_b));
                    }
                }
            }
        }

        for (i, child, new_hap_a, new_hap_b) in to_update {
            if let Some((hap_a, hap_b)) = iht_vecs[i].iht.children.get_mut(&child) {
                *hap_a = new_hap_a;
                *hap_b = new_hap_b;
            }
        }

        if !changes_made {
            break;
        }
    }
}

/// Analyzes a reference to `Vec<IhtVec>` and detects children whose characters change between consecutive segments.
/// Outputs a vector of space-separated strings in the format:
/// `last_end next_start child_id old_char new_char`
fn summarize_child_changes(iht_vecs: &Vec<IhtVec>) -> Vec<String> {
    let mut summaries = Vec::new();

    for window in iht_vecs.windows(2) {
        let previous = &window[0];
        let next = &window[1];

        for (child, &(prev_a, prev_b)) in &previous.iht.children {
            if let Some(&(next_a, next_b)) = next.iht.children.get(child) {
                // Check if either haplotype changes
                if prev_a != next_a {
                    summaries.push(format!(
                        "{} {} {} {} {}",
                        previous.bed.end, next.bed.start, child, prev_a, next_a
                    ));
                }
                if prev_b != next_b {
                    summaries.push(format!(
                        "{} {} {} {} {}",
                        previous.bed.end, next.bed.start, child, prev_b, next_b
                    ));
                }
            }
        }
    }

    summaries
}

fn optimize_and_merge_ihtvecs(
    iht_vecs: &Vec<IhtVec>,
    founders: Vec<&Individual>,
    family: &Family,
) -> Vec<IhtVec> {
    let mut optimized_vecs: Vec<IhtVec> = Vec::new();

    for mut current in iht_vecs.iter().cloned() {
        if let Some(previous) = optimized_vecs.last_mut() {
            let original_mismatches = count_mismatches(&previous.iht, &current.iht);
            let mut swapped_iht = current.iht.clone();

            for founder in &founders {
                let founder_id = founder.id();

                // Get the founder's alleles
                if let Some((founder_allele_a, founder_allele_b)) =
                    swapped_iht.founders.get(&founder_id)
                {
                    // Clone before modifying to avoid borrowing conflicts
                    let mut temp_iht = swapped_iht.clone();

                    // Swap the founder alleles across all children
                    for (child_id, (hap_a, hap_b)) in temp_iht.children.iter_mut() {
                        if *hap_a == *founder_allele_a {
                            *hap_a = *founder_allele_b;
                        } else if *hap_a == *founder_allele_b {
                            *hap_a = *founder_allele_a;
                        }

                        if *hap_b == *founder_allele_a {
                            *hap_b = *founder_allele_b;
                        } else if *hap_b == *founder_allele_b {
                            *hap_b = *founder_allele_a;
                        }
                    }

                    let swapped_mismatches = count_mismatches(&previous.iht, &temp_iht);

                    // Apply swap across children only if it reduces mismatches
                    if swapped_mismatches < original_mismatches {
                        /*
                        println!(
                            "{} {} {}/{} c:{} o:{}\no:{}\nm:{}",
                            current.bed.start,
                            founder_id,
                            founder_allele_a,
                            founder_allele_b,
                            swapped_mismatches,
                            original_mismatches,
                            current.iht.collapse_to_string(),
                            temp_iht.collapse_to_string()
                        );
                        */
                        swapped_iht = temp_iht;
                    }
                }
            }

            // Otherwise, update `current.iht`
            current.iht = swapped_iht;
        }

        optimized_vecs.push(current);
    }

    optimized_vecs
}

/// Counts the number of mismatches between two Iht structures
fn count_mismatches(iht1: &Iht, iht2: &Iht) -> usize {
    iht1.children
        .iter()
        .filter(|(child_id, (a1, b1))| {
            if let Some((a2, b2)) = iht2.children.get(*child_id) {
                (*a1 != *a2) as usize + (*b1 != *b2) as usize > 0
            } else {
                true // Mismatch if child is absent in one of the Ihts
            }
        })
        .count()
}

fn backfill_sibs(fam: &Family, iht: &Iht) -> Iht {
    let mut updated_iht = iht.clone(); // Create a mutable clone

    for (founder_id, (founder_hap_a, founder_hap_b)) in &iht.founders {
        // Get all children of this founder
        let children = fam.get_children(founder_id);

        // Process only if the founder has multiple children
        if children.len() > 1 {
            let mut identified_allele: Option<(char, usize)> = None;

            // Step 1: Identify which founder allele is present in at least one child
            for child_id in &children {
                if let Some(&(hap_a, hap_b)) = updated_iht.children.get(*child_id) {
                    if hap_a == *founder_hap_a || hap_a == *founder_hap_b {
                        identified_allele = Some((hap_a, 0));
                    } else if hap_b == *founder_hap_a || hap_b == *founder_hap_b {
                        identified_allele = Some((hap_b, 1));
                    }
                }
            }

            // Step 2: If no child has a founder allele, skip backfilling
            if identified_allele.is_none() {
                continue;
            }

            let (mut inherited_allele, allele_index) = identified_allele.unwrap();
            let mut non_inherited_allele = if inherited_allele == *founder_hap_a {
                *founder_hap_b
            } else {
                *founder_hap_a
            };

            // Step 3: Assign the other founder allele to children who do not have the inherited allele
            for child_id in &children {
                if let Some(child_hap) = updated_iht.children.get_mut(child_id.clone()) {
                    let (child_hap_a, child_hap_b) = child_hap;

                    // Check if the child already has the inherited allele
                    if *child_hap_a == inherited_allele || *child_hap_b == inherited_allele {
                        continue;
                    }

                    // Assign the missing founder allele to the correct index
                    if allele_index == 0 {
                        *child_hap_a = non_inherited_allele;
                    } else {
                        *child_hap_b = non_inherited_allele;
                    }
                }
            }

            // Step 4: Count allele frequencies AFTER backfilling
            let mut allele_counts = std::collections::HashMap::new();
            for child_id in &children {
                if let Some(&(hap_a, hap_b)) = updated_iht.children.get(*child_id) {
                    *allele_counts.entry(hap_a).or_insert(0) += 1;
                    *allele_counts.entry(hap_b).or_insert(0) += 1;
                }
            }

            let inherited_count = *allele_counts.get(&inherited_allele).unwrap_or(&0);
            let non_inherited_count = *allele_counts.get(&non_inherited_allele).unwrap_or(&0);

            // Step 5: If the inherited allele is less frequent, swap them
            if inherited_count < non_inherited_count {
                std::mem::swap(&mut inherited_allele, &mut non_inherited_allele);

                // Apply the swap across all children
                for child_id in &children {
                    if let Some((hap_a, hap_b)) = updated_iht.children.get_mut(child_id.clone()) {
                        if *hap_a == inherited_allele {
                            *hap_a = non_inherited_allele;
                        } else if *hap_a == non_inherited_allele {
                            *hap_a = inherited_allele;
                        }

                        if *hap_b == inherited_allele {
                            *hap_b = non_inherited_allele;
                        } else if *hap_b == non_inherited_allele {
                            *hap_b = inherited_allele;
                        }
                    }
                }
            }
        }
    }

    updated_iht // Return the modified Iht
}

fn find_IBD(fam: &Family, iht: &Iht) -> bool {
    for (founder_id, (founder_hap_a, founder_hap_b)) in &iht.founders {
        let children = fam.get_children(founder_id);

        // Process only if the founder has multiple children
        if children.len() > 1 {
            let mut allele_counts = std::collections::HashMap::new();

            // Step 1: Count occurrences of **only the founder alleles** in children
            for child_id in &children {
                if let Some(&(hap_a, hap_b)) = iht.children.get(*child_id) {
                    if hap_a == *founder_hap_a || hap_a == *founder_hap_b {
                        *allele_counts.entry(hap_a).or_insert(0) += 1;
                    }
                    if hap_b == *founder_hap_a || hap_b == *founder_hap_b {
                        *allele_counts.entry(hap_b).or_insert(0) += 1;
                    }
                }
            }

            // Step 2: Check if **any** founder allele appears in all children
            if allele_counts.values().any(|&count| count == children.len()) {
                return true; // IBD confirmed for this founder
            }
        }
    }

    false // No founder met the IBD condition
}

fn main() {
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

    let mut iht_vec_output_fn = args.prefix.clone();
    iht_vec_output_fn += ".iht.txt";

    let mut marker_output_fn = args.prefix.clone();
    marker_output_fn += ".markers.txt";

    let mut recomb_output_fn = args.prefix.clone();
    recomb_output_fn += ".recombinants.txt";

    let iht_path = Path::new(&iht_vec_output_fn);
    let marker_path = Path::new(&marker_output_fn);
    let recomb_path = Path::new(&recomb_output_fn);

    if iht_path.exists() || marker_path.exists() || recomb_path.exists() {
        error!("output already exists for prefix: \" {} \" ", args.prefix);
        process::exit(1); // Exit with a non-zero status code
    }

    if let Err(e) = check_individuals_in_vcf(&args.ped, &args.vcf) {
        eprintln!("Error while checking individuals in VCF: {}", e);
    }

    match is_vcf_indexed(&args.vcf) {
        Ok(true) => {}
        Ok(false) => {
            warn!("VCF file is not indexed.");
            process::exit(1); // Exit with a non-zero status code
        }
        Err(e) => {
            eprintln!("Error checking VCF index: {}", e);
            process::exit(1); // Exit with a non-zero status code for errors
        }
    }

    // Open the file with write and create options
    let mut iht_file = OpenOptions::new()
        .write(true) // Open for writing
        .create(true) // Create if it doesn't exist
        .open(iht_vec_output_fn)
        .unwrap();

    // Open the file with write and create options
    let mut marker_file = OpenOptions::new()
        .write(true) // Open for writing
        .create(true) // Create if it doesn't exist
        .open(marker_output_fn)
        .unwrap();

    // Open the file with write and create options
    let mut recomb_file = OpenOptions::new()
        .write(true) // Open for writing
        .create(true) // Create if it doesn't exist
        .open(recomb_output_fn)
        .unwrap();

    let family = Family::parse_ped_file(&args.ped).unwrap();
    let master_iht = Iht::new(family.founders(), family.offspring(), &ChromType::Autosome);

    iht_file
        .write(
            format!(
                "#chrom start end {} marker_count len\n",
                master_iht.legend()
            )
            .as_bytes(),
        )
        .unwrap();

    marker_file
        .write(format!("#chom pos founder allele matches {}\n", master_iht.legend()).as_bytes())
        .unwrap();

    recomb_file
        .write(format!("#chrom start end sample hap1 hap2\n").as_bytes())
        .unwrap();

    let mut reader: IndexedReader =
        IndexedReader::from_path(&args.vcf).expect("Failure to read VCF file.");

    let chromosomes = extract_chromosome_names(&args.vcf).unwrap();

    for c in chromosomes {
        let mut zygosity = ChromType::Autosome;

        if c.0.contains("chrX") || c.0.contains("ChrX") {
            zygosity = ChromType::ChrX;
        }

        let mut iht_info = Iht::new(family.founders(), family.offspring(), &zygosity);

        let genotype_data = parse_vcf(&mut reader, c.1, &c.0, args.qual).unwrap();

        if genotype_data.len() == 0 {
            continue;
        }

        info!("{} has {} variant records.", c.0, genotype_data.len());

        let mut pre_vector: Vec<IhtVec> = Vec::new();

        for mut gs in genotype_data {
            remove_unused_samples(&mut gs.1, &family.get_individuals_ids());

            if has_missing_alleles(&gs.1, args.depth) {
                continue;
            }

            let markers: (
                bool,
                HashMap<String, (HashSet<GenotypeAllele>, Vec<String>)>,
            ) = find_founder_unique_alleles(&gs.1, &family);
            if !markers.0 {
                continue;
            }

            let mut local_iht = Iht::new(family.founders(), family.offspring(), &zygosity);

            // going from markers to iht structure
            load_up(&family, &mut local_iht, &markers.1, &zygosity);

            // backfilling the other founder allele in kinships (multi-child families)
            let backfilled = backfill_sibs(&family, &local_iht);

            marker_file
                .write(
                    format!(
                        "{} {} {} {}\n",
                        gs.0.chrom,
                        gs.0.start,
                        marker_to_string(&markers.1),
                        backfilled.collapse_to_string()
                    )
                    .as_bytes(),
                )
                .unwrap();

            if find_IBD(&family, &backfilled) {
                warn!(
                    "IBD marker is being skipped, but will be in marker file {} {} {}",
                    gs.0.chrom,
                    gs.0.start,
                    backfilled.collapse_to_string()
                );
                continue;
            }

            iht_info.merge(&backfilled);

            pre_vector.push(IhtVec {
                bed: gs.0,
                count: 1,
                iht: backfilled.clone(),
                non_missing_counts: count_non_missing(&backfilled.founders, &backfilled.children),
            });
        }
        info!("{} has {} haplotype marker sites.", c.0, pre_vector.len());

        let flipped = optimize_and_merge_ihtvecs(
            &pre_vector,
            family.get_founders_with_multiple_children(),
            &family,
        );

        // collapsing identical marker sets
        let mut iht_vecs = collapse_identical_iht(flipped);

        // removing
        apply_threshold_to_vec(&mut iht_vecs, args.run);

        let mut filtered_iht_vec = collapse_identical_iht(iht_vecs);
        fill_missing_values(&mut filtered_iht_vec);

        let double_flipped: Vec<IhtVec> = optimize_and_merge_ihtvecs(
            &filtered_iht_vec,
            family.get_founders_with_multiple_children(),
            &family,
        );

        // collapsing identical marker sets
        let mut final_vecs = collapse_identical_iht(double_flipped);

        for i in &final_vecs {
            iht_file
                .write(
                    format!(
                        "{} {} {} {} {} {}\n",
                        i.bed.chrom,
                        i.bed.start,
                        i.bed.end,
                        i.iht.collapse_to_string(),
                        i.count,
                        i.bed.end - i.bed.start,
                    )
                    .as_bytes(),
                )
                .unwrap();
        }

        for recomb in summarize_child_changes(&final_vecs) {
            recomb_file
                .write(format!("{} {}\n", c.0, recomb).as_bytes())
                .unwrap();
        }
    }
}
