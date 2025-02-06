use clap::Parser;
use concordance::bed::BedRecord;
use concordance::iht;
use itertools::iproduct;
use itertools::Itertools;
use log::{debug, error, info, warn, LevelFilter};
use rust_htslib::htslib::fai_format_options_FAI_FASTA;
use rust_htslib::htslib::hts_name2id_f;
use std::collections::VecDeque;
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::fs::File;
use std::fs::OpenOptions;
use std::io;
use std::io::BufRead;
use std::io::Read as IoRead;
use std::io::Write;
use std::mem::swap;
use std::path::Path;
use std::process;
use std::str;

use rust_htslib::bcf::header::HeaderView;
use rust_htslib::bcf::record::{Genotype, GenotypeAllele};
use rust_htslib::bcf::Record;
use rust_htslib::bcf::{Format, Header, IndexedReader, Read, Writer};

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
#[derive(Clone, Debug)]
struct Iht {
    /// sample ID, (hapA, hapB)
    founders: HashMap<String, (char, char)>,
    children: HashMap<String, (char, char)>,
}
#[derive(Clone, Debug)]
struct IhtVec {
    bed: BedRecord,
    count: usize,
    iht: Iht,
    non_missing_counts: HashMap<String, HashMap<char, usize>>, // Store counts of each non '?' character per child
}

impl fmt::Display for Iht {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut output = String::new();

        // Sort and append founders to the output
        let mut sorted_founders: Vec<_> = self.founders.iter().collect();
        sorted_founders.sort_by_key(|(id, _)| *id);
        output.push_str("Founders:\n");
        for (id, (hap_a, hap_b)) in sorted_founders {
            output.push_str(&format!("  {} -> ({}, {})\n", id, hap_a, hap_b));
        }

        // Sort and append children to the output
        let mut sorted_children: Vec<_> = self.children.iter().collect();
        sorted_children.sort_by_key(|(id, _)| *id);
        output.push_str("Children:\n");
        for (id, (hap_a, hap_b)) in sorted_children {
            output.push_str(&format!("  {} -> ({}, {})\n", id, hap_a, hap_b));
        }

        // Write the formatted output to the formatter
        write!(f, "{}", output)
    }
}

impl Iht {
    fn new(founders: Vec<&Individual>, children: Vec<&Individual>) -> Self {
        let mut founder_defaults = ('A'..='Z').step_by(2).zip(('B'..='Z').step_by(2));

        let founders_map: HashMap<String, (char, char)> = founders
            .iter()
            .map(|id| {
                let hap_pair = founder_defaults.next().unwrap_or(('?', '?'));
                (id.id(), hap_pair)
            })
            .collect();

        let children_map: HashMap<String, (char, char)> =
            children.iter().map(|id| (id.id(), ('?', '?'))).collect();

        Self {
            founders: founders_map,
            children: children_map,
        }
    }

    /// Checks both founders and children for a key and returns the value if found.
    fn get_alleles(&self, key: &str) -> Option<(char, char)> {
        self.founders
            .get(key)
            .cloned()
            .or_else(|| self.children.get(key).cloned())
    }

    fn generate_combinations(&self, family: &Family) -> Vec<Iht> {
        let mut results: Vec<Iht> = Vec::new();

        fn generate_tuples(n: usize) -> Vec<Vec<(u8, u8)>> {
            if n == 0 {
                return vec![]; // Return an empty vector if no samples
            }

            // Define the base tuples
            let base = vec![(0_u8, 1_u8), (1_u8, 0_u8)];

            // Create `n` iterators of the base tuples `(0, 1)` and `(1, 0)`
            let iterators = std::iter::repeat(base.into_iter()).take(n);

            // Generate the Cartesian product of `n` samples
            iterators
                .multi_cartesian_product() // Cartesian product across all dimensions
                .collect::<Vec<Vec<(u8, u8)>>>()
        }

        fn generate_founder_permutations(
            founders: &HashMap<String, (char, char)>,
        ) -> Vec<HashMap<String, (char, char)>> {
            // Generate all possible allele swaps for each founder
            let mut permutations = vec![founders.clone()];

            for (id, (allele_a, allele_b)) in founders {
                let mut swapped = founders.clone();
                swapped.insert(id.clone(), (*allele_b, *allele_a)); // Swap alleles
                permutations.push(swapped);
            }

            permutations
        }

        // Generate all combinations based on the number of children
        let combs = generate_tuples(self.children.len());

        // Generate all founder allele permutations
        let founder_permutations = generate_founder_permutations(&self.founders);

        // Use `get_individual_depths` to iterate over individuals in increasing depth
        let sorted_individuals = family.get_individual_depths();

        // Iterate over all founder permutations
        for permuted_founders in founder_permutations {
            // Iterate over all combinations of allele patterns
            for com in &combs {
                let mut current_iht = self.clone();
                current_iht.founders = permuted_founders.clone(); // Update founders with the permuted alleles
                let mut child_counter = 0;

                // Iterate over individuals sorted by depth
                for (sample_id, _) in &sorted_individuals {
                    // Skip founders and non-children
                    if !self.children.contains_key(sample_id) {
                        continue;
                    }

                    // Get the child pattern for the current individual
                    if child_counter >= com.len() {
                        warn!(
                            "Child counter {} exceeds allele pattern length {} for sample: {}",
                            child_counter,
                            com.len(),
                            sample_id
                        );
                        break; // Stop processing this combination
                    }

                    if let Some((c1, c2)) = com.get(child_counter) {
                        child_counter += 1;

                        // Get the parents of the individual
                        if let Some((father, mother)) = family.get_parents(sample_id) {
                            // Get alleles for the parents
                            if let (Some((f1, f2)), Some((m1, m2))) = (
                                current_iht.get_alleles(&father),
                                current_iht.get_alleles(&mother),
                            ) {
                                // Determine child's alleles based on the pattern
                                let father_allele = if *c1 == 0 { f1 } else { f2 };
                                let mother_allele = if *c2 == 0 { m1 } else { m2 };

                                // Assign the child's alleles in `current_iht`
                                current_iht
                                    .children
                                    .insert(sample_id.clone(), (father_allele, mother_allele));
                            } else {
                                warn!(
                                    "Alleles not found for parents of child: {} (Father: {}, Mother: {})",
                                    sample_id, father, mother
                                );
                                current_iht.children.insert(sample_id.clone(), ('?', '?'));
                                // Fallback
                            }
                        } else {
                            warn!(
                                "Parents not found for child: {} (but listed in `self.children`)",
                                sample_id
                            );
                            current_iht.children.insert(sample_id.clone(), ('?', '?'));
                            // Fallback
                        }
                    } else {
                        warn!(
                            "No allele pattern found for child: {} at position {}",
                            sample_id, child_counter
                        );
                        current_iht.children.insert(sample_id.clone(), ('?', '?'));
                        // Fallback
                    }
                }

                // Add the updated `current_iht` to the results
                results.push(current_iht);
            }
        }

        results
    }

    fn legend(&self) -> String {
        let mut result = String::new();

        // Sort and add founders to the string
        let mut sorted_founders: Vec<_> = self.founders.iter().map(|n| n.0.clone()).collect();
        sorted_founders.sort();

        result += &sorted_founders.join(" ");
        result += " ";

        // Sort and add children to the string
        let mut sorted_children: Vec<_> = self.children.iter().map(|n| n.0.clone()).collect();
        sorted_children.sort();
        result += &sorted_children.join(" ");

        result
    }

    /// Collapse the Iht structure into a sorted string representation.
    ///
    /// The resulting string concatenates the two characters for each individual,
    /// sorted by sample ID (founders first, then children).
    fn collapse_to_string(&self) -> String {
        // Sort and add founders to the string
        let mut sorted_founders: Vec<_> = self.founders.iter().collect();
        sorted_founders.sort_by_key(|(id, _)| *id);
        let founders_str: String = sorted_founders
            .iter()
            .map(|(_, (hap_a, hap_b))| format!("{}/{} ", hap_a, hap_b))
            .collect::<Vec<_>>()
            .join(" ");

        // Sort and add children to the string
        let mut sorted_children: Vec<_> = self.children.iter().collect();
        sorted_children.sort_by_key(|(id, _)| *id);
        let children_str: String = sorted_children
            .iter()
            .map(|(_, (hap_a, hap_b))| format!("{}|{} ", hap_a, hap_b))
            .collect::<Vec<_>>()
            .join(" ");

        format!("{} {}", founders_str, children_str)
    }

    /// Merges another `Iht` object into this one, following simplified merge rules.
    ///
    /// # Merge Rules:
    /// - Only children are merged (founders remain unchanged).
    /// - `"?"` is a placeholder that can be replaced by any non-"?" character.
    /// - If both alleles are different and non-"?", the second one overwrites the first.
    /// - **If `hap_a == hap_b` after merging, the function panics immediately, printing both original structures**.
    ///
    /// # Panics:
    /// - If after merging, a child has the same value for both haplotypes (`hap_a == hap_b`).
    fn merge(&mut self, other: &Iht) {
        let original_self = self.clone(); // Backup of the original structure before merging

        for (child_id, (hap_a_other, hap_b_other)) in &other.children {
            if let Some((hap_a_self, hap_b_self)) = self.children.get_mut(child_id) {
                // Apply simplified merge rules
                *hap_a_self = Self::resolve_merge(*hap_a_self, *hap_a_other);
                *hap_b_self = Self::resolve_merge(*hap_b_self, *hap_b_other);

                // Check if hap_a == hap_b after merge (ERROR condition)
                if *hap_a_self == *hap_b_self && *hap_a_self != '?' {
                    eprintln!(
                        "ERROR: Merge resulted in identical alleles for child `{}`",
                        child_id
                    );
                    eprintln!("Original Iht (before merge attempt): {}", original_self);
                    eprintln!("Other Iht (being merged): {}", other);
                    panic!(
                        "Merge resulted in invalid Iht where hap_a == hap_b for child `{}`",
                        child_id
                    );
                }
            } else {
                // If the child is not present in self, add it from other
                self.children
                    .insert(child_id.clone(), (*hap_a_other, *hap_b_other));
            }
        }
    }

    /// Resolves the merge of two alleles according to the new rules.
    ///
    /// # Rules:
    /// - `"?"` is replaced by any non-"?" character.
    /// - If both alleles are different and non-"?", the second allele overwrites the first.
    ///
    /// # Returns:
    /// - The merged allele.
    fn resolve_merge(allele_self: char, allele_other: char) -> char {
        match (allele_self, allele_other) {
            ('?', x) => x, // Replace "?" with actual allele
            (x, '?') => x, // Replace "?" with actual allele
            (_, y) => y,   // Overwrite the first allele with the second
        }
    }
}

fn check_individuals_in_vcf(ped_path: &str, vcf_path: &str) -> io::Result<()> {
    let family = Family::parse_ped_file(ped_path)?;

    let vcf_samples = get_samples(vcf_path).unwrap();

    for id in family.get_individuals_ids() {
        if !vcf_samples.contains(&id) {
            warn!(
                "Individual ID {} from PED file is not found in the VCF header.",
                id
            );
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
    let mut bcf = IndexedReader::from_path(vcf_path).expect("Error opening vcf file.");
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
    let mut header = reader.header();
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
) {
    for f in markers {
        let founder_alleles = iht.get_alleles(f.0).unwrap();
        let founder_info = fam.get_individual(f.0).unwrap();
        for o in &f.1 .1 {
            let child_alleles = iht.children.get(o).unwrap();
            // paternal goes left

            // this bit of logic deals with fact that founders alleles are passed across generations and we need to flip past G2, based on the parent.
            let child = fam.get_individual(&o).unwrap();
            let father = fam.get_individual(&child.get_father_id().unwrap()).unwrap();

            if f.1 .1.contains(&father.id()) || *f.0 == father.id() {
                iht.children
                    .insert((*o).clone(), (founder_alleles.0, child_alleles.1));
            } else {
                // maternal goes right
                iht.children
                    .insert((*o).clone(), (child_alleles.1, founder_alleles.0));
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

impl PartialEq for Iht {
    fn eq(&self, other: &Self) -> bool {
        if self.founders.len() != other.founders.len()
            || self.children.len() != other.children.len()
        {
            return false;
        }

        self.founders.iter().all(|(key, &(a1, a2))| {
            other.founders.get(key).map_or(false, |&(b1, b2)| {
                (a1 == b1 || a1 == '?' || b1 == '?') && (a2 == b2 || a2 == '?' || b2 == '?')
            })
        }) && self.children.iter().all(|(key, &(a1, a2))| {
            other.children.get(key).map_or(false, |&(b1, b2)| {
                (a1 == b1 || a1 == '?' || b1 == '?') && (a2 == b2 || a2 == '?' || b2 == '?')
            })
        })
    }
}

fn collapse_identical_iht(data: Vec<IhtVec>) -> Vec<IhtVec> {
    let mut collapsed = Vec::new();
    let mut iter = data.into_iter().peekable();

    while let Some(mut current) = iter.next() {
        let mut start = current.bed.start;
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

fn get_sample_depths(record: &Record, samples: &Vec<String>) -> Option<HashMap<String, i32>> {
    let dp_values = record.format(b"DP").integer().ok()?; // Use `ok()?` to handle potential errors

    let mut depths = HashMap::new();
    for (i, dp_sample) in dp_values.iter().enumerate() {
        if let Some(&depth) = dp_sample.get(0) {
            if let Some(sample_name) = samples.get(i) {
                depths.insert(sample_name.clone(), depth); // Clone the sample name to store as owned String
            }
        }
    }

    Some(depths) // Wrap in `Some()` to match the function signature
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
    iht_vec_output_fn += ".iht";

    let mut marker_output_fn = args.prefix.clone();
    marker_output_fn += ".markers";

    let iht_path = Path::new(&iht_vec_output_fn);
    let marker_path = Path::new(&marker_output_fn);

    if iht_path.exists() || marker_path.exists() {
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

    let family = Family::parse_ped_file(&args.ped).unwrap();
    let mut iht_info = Iht::new(family.founders(), family.offspring());

    iht_file
        .write(format!("#chrom start end {} marker_count len\n", iht_info.legend()).as_bytes())
        .unwrap();

    marker_file
        .write(format!("#chom pos founder allele matches {}\n", iht_info.legend()).as_bytes())
        .unwrap();

    let mut reader: IndexedReader =
        IndexedReader::from_path(&args.vcf).expect("Failure to read VCF file.");

    let chromosomes = extract_chromosome_names(&args.vcf).unwrap();

    for c in chromosomes {
        let mut genotype_data = parse_vcf(&mut reader, c.1, &c.0, args.qual).unwrap();

        if genotype_data.len() == 0 {
            continue;
        }

        info!(
            "{} has {} variant records (sites)",
            c.0,
            genotype_data.len()
        );

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

            let mut local_iht = Iht::new(family.founders(), family.offspring());

            load_up(&family, &mut local_iht, &markers.1);

            marker_file
                .write(
                    format!(
                        "{} {} {} {}\n",
                        gs.0.chrom,
                        gs.0.start,
                        marker_to_string(&markers.1),
                        local_iht.collapse_to_string()
                    )
                    .as_bytes(),
                )
                .unwrap();

            iht_info.merge(&local_iht);

            pre_vector.push(IhtVec {
                bed: gs.0,
                count: 1,
                iht: local_iht.clone(),
                non_missing_counts: count_non_missing(&local_iht.founders, &local_iht.children),
            });
        }
        info!("There are {} haplotype marker (sites).", pre_vector.len());

        let mut iht_vecs = collapse_identical_iht(pre_vector);

        apply_threshold_to_vec(&mut iht_vecs, args.run);

        let mut filtered_iht_vec = collapse_identical_iht(iht_vecs);
        fill_missing_values(&mut filtered_iht_vec);

        for i in filtered_iht_vec {
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
    }
}
