use clap::Parser;
use concordance::bed::BedRecord;
use concordance::iht::ChromType;
use concordance::iht::Iht;
use concordance::iht::IhtVec;
use concordance::utils::*;

use itertools::Itertools;
use log::{debug, error, info, warn, LevelFilter};

use std::collections::{HashMap, HashSet};

use std::fs::OpenOptions;
use std::io;
use std::io::Write;

use std::collections::VecDeque;
use std::path::Path;
use std::process;
use std::str;
use std::usize;

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
            let unphased = unphase_genotype(&alleles);

            //println!("before: {:?} after {:?}", alleles, unphased);
            record_map.insert(sample.clone(), (*depths.get(sample).unwrap(), unphased));
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

/// Checks if any sample has missing alleles, a sequencing depth of zero,
/// or an abnormally high depth (>2 std deviations above sample mean).
///
/// # Arguments
/// - `genotype_data`: A `HashMap` mapping sample IDs to their depth and genotype alleles.
/// - `depth_stats`: A `HashMap` from `extract_depth_statistics` containing sample mean and std dev.
/// - `min_depth`: A minimum depth threshold (e.g., 10).
///
/// # Returns
/// - `true` if any sample has a depth below `min_depth`, is missing alleles,
///   or has a depth > 2 std deviations above the sample's mean or below.
/// - `false` otherwise.
fn depth_filters(
    genotype_data: &HashMap<String, (i32, Vec<GenotypeAllele>)>,
    depth_stats: &HashMap<String, (f64, f64)>,
    min_depth: i32,
) -> bool {
    // First pass: Check for missing alleles or depth below min_depth
    for (_sample, (depth, allele_vec)) in genotype_data.iter() {
        if *depth < min_depth || allele_vec.iter().any(|allele| allele.index().is_none()) {
            return true;
        }
    }

    // Second pass: Check for extreme depth values (more than 2 standard deviations from mean)
    for (sample, (depth, _)) in genotype_data.iter() {
        if let Some((mean, std_dev)) = depth_stats.get(sample) {
            let lower_bound = mean - (1.0 * std_dev);
            let upper_bound = mean + (1.0 * std_dev);
            if (*depth as f64) < lower_bound || (*depth as f64) > upper_bound {
                return true;
            }
        }
    }

    false
}

fn unique_allele(
    individual1: &Vec<GenotypeAllele>,
    individual2: &Vec<GenotypeAllele>,
    zygosity: &ChromType,
    sex: u8,
) -> Option<GenotypeAllele> {
    // if the chromosome is X and the individual 1 is female, she needs to be het (valid marker)
    if (*zygosity == ChromType::ChrX && sex == 2) || (*zygosity != ChromType::ChrX) {
        if individual1.get(0).unwrap() == individual1.get(1).unwrap() {
            return None;
        }
    }

    // if the chromosome is X and the individual is 1 male, he needs to be hom (deepVariant)
    if *zygosity == ChromType::ChrX && sex == 1 {
        if individual1.get(0).unwrap() != individual1.get(1).unwrap() {
            return None;
        }
    }

    let set1: HashSet<GenotypeAllele> = individual1.iter().cloned().collect();
    let set2: HashSet<GenotypeAllele> = individual2.iter().cloned().collect();

    // Find the allele in individual1 that is not in individual2
    return set1.difference(&set2).next().cloned();
}

fn has_allele(a: &GenotypeAllele, geno: &Vec<GenotypeAllele>) -> bool {
    geno.iter().any(|b| b == a)
}

fn get_iht_markers(sample: &String, local_iht: &mut Iht) -> (char, char) {
    if local_iht.children.contains_key(sample) {
        return *local_iht.children.get(sample).unwrap();
    }
    if local_iht.founders.contains_key(sample) {
        return *local_iht.founders.get(sample).unwrap();
    }

    ('?', '?')
}

fn find_valid_char(chars: (char, char)) -> Option<char> {
    if chars.0 != '?' && chars.0 != '.' {
        Some(chars.0)
    } else if chars.1 != '?' && chars.1 != '.' {
        Some(chars.1)
    } else {
        None
    }
}

fn track_alleles_through_pedigree(
    genotype_data: &HashMap<String, (i32, Vec<GenotypeAllele>)>,
    family: &Family,
    local_iht: &mut Iht,
    zygosity: &ChromType,
) -> (
    bool,
    HashMap<String, (HashSet<GenotypeAllele>, Vec<String>)>,
) {
    let mut inherited_marker_allele: HashMap<String, GenotypeAllele> = HashMap::new();

    // walk down the family from founders to offspring
    for ind_id in family.get_individual_depths() {
        let ind = family.get_individual(&ind_id.0).unwrap();

        let (_, ind_alleles) = genotype_data.get(&ind.id()).unwrap();
        // individual has a spouse
        if let Some(spouse_id) = family.find_spouse(&ind.id()) {
            let (_, spouse_alleles) = genotype_data.get(&spouse_id).unwrap();
            // the individual has a unique trackable marker
            if let Some(marker) = unique_allele(
                ind_alleles,
                spouse_alleles,
                zygosity,
                ind.get_sex().unwrap(),
            ) {
                if let Some(iht_marker) = inherited_marker_allele.get(&ind_id.0) {
                    // if the unique allele isn't the marker we are tracking then we need to skip it.
                    if marker != *iht_marker {
                        continue;
                    }
                }

                let flabels = get_iht_markers(&ind_id.0, local_iht);

                if flabels.0 == flabels.1 && flabels.0 == '?' {
                    continue;
                }
                let flabel = find_valid_char(flabels).unwrap();

                // looping over children
                for c in family.get_children(&ind.id()) {
                    let child_fam = family.get_individual(c).unwrap();
                    if has_allele(&marker, &genotype_data.get(c).unwrap().1) {
                        #[allow(unused_assignments)]
                        let child = local_iht.children.get_mut(c).unwrap();

                        let parent_sex = ind.get_sex().unwrap();
                        let child_sex = child_fam.get_sex().unwrap();

                        // If the child is a male it cannot get X from dad, so we continue
                        if child_sex == 1 && parent_sex == 1 && *zygosity == ChromType::ChrX {
                            continue;
                        }
                        // This includes dad passing X to daughter
                        if ind.get_sex().unwrap() == 1 {
                            // parent is male so we are passing to the first position
                            child.0 = flabel;
                        } else {
                            // parent is female so we are passing to the first position
                            child.1 = flabel;
                        }

                        inherited_marker_allele.insert(c.clone(), marker);
                    }
                }
            }
        }
    }

    // After processing markers, update children genotypes for ChrX:
    // If the ChromType is ChrX and the child is male, set the first genotype entry to "."
    if *zygosity == ChromType::ChrX {
        for (child_id, genotype) in local_iht.children.iter_mut() {
            let child = family.get_individual(child_id).unwrap();
            if child.get_sex().unwrap() == 1 {
                genotype.0 = '.';
            }
        }
    }

    for c in &local_iht.children {
        if (c.1 .0 != '?' && c.1 .0 != '.') || (c.1 .1 != '?' && c.1 .1 != '.') {
            return (true, HashMap::new());
        }
    }

    return (false, HashMap::new());
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
                debug!("merging iht blocks");
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
fn fill_missing_values_by_neighbor(iht_vecs: &mut Vec<IhtVec>) {
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

/// Fills in '?' values when matching entries exist before and after, iterating until no more changes occur.
fn fill_missing_values(iht_vecs: &mut Vec<IhtVec>, pre_filter: Vec<IhtVec>) {
    for iht_vec in iht_vecs.iter_mut() {
        for (child, (hap_a, hap_b)) in iht_vec.iht.children.iter_mut() {
            let mut consensus_a = None;
            let mut consensus_b = None;

            if *hap_a == '?' || *hap_b == '?' {
                let mut marker_counts_a: HashMap<char, usize> = HashMap::new();
                let mut marker_counts_b: HashMap<char, usize> = HashMap::new();

                for pre in pre_filter.iter() {
                    if pre.bed.start >= iht_vec.bed.start && pre.bed.end <= iht_vec.bed.end {
                        if let Some((pre_a, pre_b)) = pre.iht.children.get(child) {
                            if *pre_a != '?' {
                                *marker_counts_a.entry(*pre_a).or_insert(0) += 1;
                            }
                            if *pre_b != '?' {
                                *marker_counts_b.entry(*pre_b).or_insert(0) += 1;
                            }
                        }
                    }
                }

                if !marker_counts_a.is_empty() {
                    consensus_a = marker_counts_a
                        .iter()
                        .max_by_key(|entry| entry.1)
                        .map(|(k, _)| *k);
                }
                if !marker_counts_b.is_empty() {
                    consensus_b = marker_counts_b
                        .iter()
                        .max_by_key(|entry| entry.1)
                        .map(|(k, _)| *k);
                }

                if *hap_a == '?' {
                    if let Some(consensus) = consensus_a {
                        *hap_a = consensus;
                    }
                }

                if *hap_b == '?' {
                    if let Some(consensus) = consensus_b {
                        *hap_b = consensus;
                    }
                }
            }
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

fn perform_flips_in_place(
    iht_vecs: &mut Vec<IhtVec>,
    founders: Vec<&Individual>,
    family: &Family,
    zygosity: &ChromType,
) {
    if *zygosity == ChromType::ChrX {
        return;
    }

    let flippable_indices: Vec<usize> = iht_vecs
        .iter()
        .enumerate()
        .filter_map(|(idx, iht_vec)| {
            if !iht_vec.iht.get_flipable_alleles(family).is_empty() {
                Some(idx)
            } else {
                None
            }
        })
        .collect();

    for i in 0..iht_vecs.len() {
        let current_flipable = iht_vecs[i].iht.get_flipable_alleles(family);
        if current_flipable.is_empty() {
            continue;
        }

        // Find the most recent previous index with flippable alleles
        let prev_idx = flippable_indices.iter().rev().find(|&&idx| idx < i);

        if let Some(&prev_i) = prev_idx {
            // Using `get_mut` to safely obtain a mutable reference without double borrowing
            let (previous, current) = iht_vecs.split_at_mut(prev_i + 1);
            let previous = &mut previous[prev_i];
            let current = &mut current[i - prev_i - 1];

            let prev_flipable = previous.iht.get_flipable_alleles(family);

            if prev_flipable == current_flipable
                || prev_flipable.is_superset(&current_flipable)
                || current_flipable.is_superset(&prev_flipable)
            {
                let original_mismatches = count_mismatches(&previous.iht, &current.iht);
                let mut swapped_iht = current.iht.clone();

                for founder in &founders {
                    if original_mismatches == 0 {
                        break;
                    }

                    let founder_id = founder.id();

                    // Get the founder's alleles
                    if let Some((founder_allele_a, founder_allele_b)) =
                        swapped_iht.founders.get(&founder_id)
                    {
                        let mut temp_iht = swapped_iht.clone();

                        // Swap the founder alleles across all children
                        for (_, (hap_a, hap_b)) in temp_iht.children.iter_mut() {
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
                            swapped_iht = temp_iht;
                        }
                    }
                }
                current.iht = swapped_iht;
            }
        }
    }
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

fn backfill_sibs(fam: &Family, iht: &Iht, zygosity: &ChromType) -> Iht {
    let mut updated_iht = iht.clone(); // Create a mutable clone

    for (founder_id, (founder_hap_a, founder_hap_b)) in &iht.founders {
        let founder = fam.get_individual(&founder_id).unwrap();

        if *zygosity == ChromType::ChrX && founder.get_sex().unwrap() == 1 {
            continue;
        }

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

            let mut iht_allele_children: HashSet<String> = HashSet::new();
            let mut non_iht_allele_children: HashSet<String> = HashSet::new();

            // Step 3: Assign the other founder allele to children who do not have the inherited allele
            for child_id in &children {
                if let Some(child_hap) = updated_iht.children.get_mut(*child_id) {
                    let (child_hap_a, child_hap_b) = child_hap;

                    // Check if the child already has the inherited allele
                    if *child_hap_a == inherited_allele || *child_hap_b == inherited_allele {
                        iht_allele_children.insert((*child_id).clone());
                        continue;
                    }

                    non_iht_allele_children.insert((*child_id).clone());

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
            if inherited_count < non_inherited_count
                || (inherited_count == non_inherited_count)
                    && iht_allele_children.iter().min().unwrap()
                        > non_iht_allele_children.iter().min().unwrap()
            {
                std::mem::swap(&mut inherited_allele, &mut non_inherited_allele);

                // Apply the swap across all children
                for child_id in &children {
                    if let Some((hap_a, hap_b)) = updated_iht.children.get_mut(*child_id) {
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
    /*
            println!("legend: {}", iht.legend());
            println!("before: {}", iht.collapse_to_string());
            println!("after:  {}\n", updated_iht.collapse_to_string());
    */
    updated_iht // Return the modified Iht
}

fn collect_alleles_with_positions(
    iht_vecs: &Vec<IhtVec>,
    sample: &str,
    index: usize,
) -> VecDeque<(i64, char)> {
    let mut allele_positions = VecDeque::new();

    for iht_vec in iht_vecs {
        if let Some((hap_a, hap_b)) = iht_vec.iht.children.get(sample) {
            let allele = if index == 0 { *hap_a } else { *hap_b };
            if allele != '?' {
                allele_positions.push_back((iht_vec.bed.start, allele));
            }
        }
    }

    allele_positions
}

fn count_matching_neighbors(
    allele_positions: &VecDeque<(i64, char)>,
    min_run: usize,
) -> Vec<(i64, char, usize, usize)> {
    let mut results = Vec::new();

    for (i, &(pos, allele)) in allele_positions.iter().enumerate() {
        let mut count_before = 1;
        let mut count_after = 1;

        // Count matches before
        for j in (0..i).rev() {
            if allele_positions[j].1 == allele {
                count_before += 1;
            } else {
                break;
            }
        }

        // Count matches after
        for j in (i + 1)..allele_positions.len() {
            if allele_positions[j].1 == allele {
                count_after += 1;
            } else {
                break;
            }
        }
        if count_before < min_run && count_after < min_run {
            results.push((pos, allele, count_before, count_after));
        }
    }

    results
}

fn mask_child_alleles(
    sites: &HashSet<i64>,
    child_id: &str,
    allele_index: usize,
    iht_vecs: &mut Vec<IhtVec>,
) {
    for iht_vec in iht_vecs.iter_mut() {
        if sites.contains(&iht_vec.bed.start) {
            if let Some((hap_a, hap_b)) = iht_vec.iht.children.get_mut(child_id) {
                if allele_index == 0 {
                    *hap_a = '?';
                } else if allele_index == 1 {
                    *hap_b = '?';
                }
            }
        }
    }
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
                "#chrom start end {} marker_count len markers\n",
                master_iht.legend()
            )
            .as_bytes(),
        )
        .unwrap();

    marker_file
        .write(format!("#chom pos {} info\n", master_iht.legend()).as_bytes())
        .unwrap();

    recomb_file
        .write(format!("#chrom start end sample hap1 hap2\n").as_bytes())
        .unwrap();

    let mut reader: IndexedReader =
        IndexedReader::from_path(&args.vcf).expect("Failure to read VCF file.");

    let chromosomes = extract_chromosome_names(&args.vcf).unwrap();

    for c in chromosomes {
        let depth_info = extract_depth_statistics(&mut reader, c.1).unwrap();
        let mut zygosity = ChromType::Autosome;

        if c.0.contains("chrX") || c.0.contains("ChrX") {
            zygosity = ChromType::ChrX;
        }

        let genotype_data = parse_vcf(&mut reader, c.1, &c.0, args.qual).unwrap();

        if genotype_data.len() == 0 {
            continue;
        }

        info!("{} has {} variant records.", c.0, genotype_data.len());

        let mut pre_vector: Vec<IhtVec> = Vec::new();

        let mut marker_info = HashMap::new();

        for mut gs in genotype_data {
            remove_unused_samples(&mut gs.1, &family.get_individuals_ids());

            if depth_filters(&gs.1, &depth_info, args.depth) {
                continue;
            }

            let mut local_iht = Iht::new(family.founders(), family.offspring(), &zygosity);

            let markers: (
                bool,
                HashMap<String, (HashSet<GenotypeAllele>, Vec<String>)>,
            ) = track_alleles_through_pedigree(&gs.1, &family, &mut local_iht, &zygosity);
            if !markers.0 {
                continue;
            }

            // backfilling the other founder allele in kinships (multi-child families)
            let backfilled = backfill_sibs(&family, &local_iht, &zygosity);

            marker_info.insert(gs.0.start, marker_to_string(&markers.1));

            pre_vector.push(IhtVec {
                bed: gs.0,
                count: 1,
                iht: backfilled.clone(),
                non_missing_counts: count_non_missing(&backfilled.founders, &backfilled.children),
            });
        }
        info!("{} has {} haplotype marker sites.", c.0, pre_vector.len());

        perform_flips_in_place(
            &mut pre_vector,
            family.get_founders_with_multiple_children(),
            &family,
            &zygosity,
        );

        for m in &pre_vector {
            let child_alleles = m.iht.get_non_missing_child_alleles();
            let mut assignment: Vec<String> = Vec::new();

            for a in child_alleles {
                let children = m.iht.get_children_by_allele(a);
                assignment.push(format!("{}:{}", a, children.join(",")));
            }

            marker_file
                .write(
                    format!(
                        "{} {} {} {}\n",
                        m.bed.chrom,
                        m.bed.start,
                        m.iht.collapse_to_string(),
                        assignment.join(";")
                    )
                    .as_bytes(),
                )
                .unwrap();
        }

        let stable_iht = pre_vector.clone();

        info!("Finding short marker runs to mask.");
        for d in family.offspring() {
            for i in [0 as usize, 1 as usize].iter() {
                let alleles = collect_alleles_with_positions(&pre_vector, &d.id(), *i);
                let neighbor_counts: Vec<(i64, char, usize, usize)> =
                    count_matching_neighbors(&alleles, args.run);
                let mut to_mask: HashSet<i64> = HashSet::new();
                let mut mask_print: Vec<i64> = Vec::new();

                for (key, _, _, _) in neighbor_counts {
                    to_mask.insert(key);
                    mask_print.push(key);
                }

                debug!(
                    "masking {} family member haplotype {} positions {}",
                    d.id(),
                    i,
                    mask_print.iter().map(ToString::to_string).join(",")
                );
                mask_child_alleles(&to_mask, &d.id(), *i, &mut pre_vector);
            }
        }

        let mut iht_vecs = collapse_identical_iht(pre_vector);

        perform_flips_in_place(
            &mut iht_vecs,
            family.get_founders_with_multiple_children(),
            &family,
            &zygosity,
        );

        fill_missing_values(&mut iht_vecs, stable_iht);
        fill_missing_values_by_neighbor(&mut iht_vecs);

        perform_flips_in_place(
            &mut iht_vecs,
            family.get_founders_with_multiple_children(),
            &family,
            &zygosity,
        );

        for i in &iht_vecs {
            iht_file
                .write(
                    format!(
                        "{} {} {} {} {} {} {}\n",
                        i.bed.chrom,
                        i.bed.start,
                        i.bed.end,
                        i.iht.collapse_to_string(),
                        i.count,
                        i.bed.end - i.bed.start,
                        i.iht.get_non_missing_child_alleles().iter().join(","),
                    )
                    .as_bytes(),
                )
                .unwrap();
        }

        for recomb in summarize_child_changes(&iht_vecs) {
            recomb_file
                .write(format!("{} {}\n", c.0, recomb).as_bytes())
                .unwrap();
        }
    }
}
