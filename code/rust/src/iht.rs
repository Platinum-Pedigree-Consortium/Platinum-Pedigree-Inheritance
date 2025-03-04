use crate::bed::BedRecord;
use crate::ped::Family;
use crate::ped::Individual;
use core::fmt;
use csv::ReaderBuilder;
use itertools::Itertools;
use log::warn;
use rust_htslib::bcf::record::GenotypeAllele;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

#[derive(PartialEq)]
pub enum ChromType {
    ChrX,
    ChrY,
    ChrM,
    Autosome,
}

pub struct InheritanceBlock {
    pub chrom: String,
    pub start: i32,
    pub end: i32,
    pub passing_count: i32,
    pub failing_count: i32,
    pub samples: Vec<String>,
    pub sample_lookups: HashMap<String, usize>,
    pub parental_hap: Vec<String>,
    pub patterns: HashMap<String, Vec<[i32; 2]>>,
    pub inherited_haps: HashSet<char>,
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

pub fn get_iht_block<'a>(
    ihts: &'a mut Vec<InheritanceBlock>,
    chr: &str,
    pos: i32,
    current: &mut usize,
) -> Option<&'a mut InheritanceBlock> {
    for i in *current..ihts.len() {
        if (ihts[i].chrom == chr) && (pos >= ihts[i].start) && (pos <= ihts[i].end) {
            return Some(&mut ihts[i]);
        }
        *current += 1;
    }
    return None;
}

pub fn parse_inht(inht_fn: String) -> Vec<InheritanceBlock> {
    use std::fs::File;
    use std::io::Seek;

    let mut inht_fp = File::open(&inht_fn).expect("Error reading inheritance CSV file.");
    let inht_fp_gz = flate2::read::GzDecoder::new(&mut inht_fp);
    let inht_fp: Box<dyn std::io::Read> = match inht_fp_gz.header() {
        Some(_) => Box::new(inht_fp_gz),
        None => {
            inht_fp.rewind().unwrap();
            Box::new(inht_fp)
        }
    };
    let mut reader = ReaderBuilder::new().from_reader(inht_fp);

    let mut inht_info = Vec::new();
    let header = reader
        .headers()
        .expect("Error reading inheritance CSV header")
        .clone();

    for record in reader.records() {
        let mut ihtblock = InheritanceBlock {
            chrom: record.as_ref().unwrap()[0].to_string().clone(),
            start: record.as_ref().unwrap()[1].parse::<i32>().unwrap().clone(),
            end: record.as_ref().unwrap()[2].parse::<i32>().unwrap().clone(),
            passing_count: 0,
            failing_count: 0,
            samples: Vec::new(),
            sample_lookups: HashMap::new(),
            parental_hap: Vec::new(),
            patterns: HashMap::new(),
            inherited_haps: HashSet::new(),
        };

        let mut one_parent = false;

        let mut sidx: usize = 0;
        for i in 3..header.len() {
            // Some inheritance vectors do not contain both parental marker, for example you might see `A` rather than `AB`.
            // We want to skip these sites as the expected genotypes are unknown.
            if record.as_ref().unwrap()[i].to_string().len() == 1 {
                one_parent = true;
            }
            // putting the sample names in the header into the inheritance block for easy lookup
            ihtblock
                .sample_lookups
                .insert((&header[i]).to_string(), sidx);
            ihtblock.samples.push((&header[i]).to_string());
            sidx += 1;
            // geno
            ihtblock
                .parental_hap
                .push(record.as_ref().unwrap()[i].to_string());
        }
        // counting up haplotypes seen in children
        for i in 0..ihtblock.parental_hap.len() - 2 {
            let geno = ihtblock.parental_hap.get(i).unwrap();
            ihtblock.inherited_haps.insert(geno.as_bytes()[0].into());
            ihtblock.inherited_haps.insert(geno.as_bytes()[1].into());
        }

        if one_parent {
            println!("Warning skipping block missing both parents {}", ihtblock);
            continue;
        }
        inht_info.push(ihtblock);
    }
    return inht_info;
}

#[derive(Clone, Debug)]
pub struct Iht {
    /// sample ID, (hapA, hapB)
    pub founders: HashMap<String, (char, char)>,
    pub children: HashMap<String, (char, char)>,
}
#[derive(Clone, Debug)]
pub struct IhtVec {
    pub bed: BedRecord,
    pub count: usize,
    pub iht: Iht,
    pub non_missing_counts: HashMap<String, HashMap<char, usize>>, // Store counts of each non '?' character per child
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
    pub fn new(founders: Vec<&Individual>, children: Vec<&Individual>, hemi: &ChromType) -> Self {
        let mut founder_defaults = ('A'..='Z').step_by(2).zip(('B'..='Z').step_by(2));

        let founders_map: HashMap<String, (char, char)> = founders
            .iter()
            .map(|id| {
                // mitochondria and Y are always hemi
                let mut hap_pair = founder_defaults.next().unwrap_or(('?', '?'));
                if *hemi == ChromType::ChrY || *hemi == ChromType::ChrM {
                    hap_pair = (hap_pair.0, '.');
                    // male chr X
                } else if *hemi == ChromType::ChrX && id.get_sex().unwrap() == 1 {
                    hap_pair = ('.', hap_pair.0);
                }

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
    pub fn get_alleles(&self, key: &str) -> Option<(char, char)> {
        self.founders
            .get(key)
            .cloned()
            .or_else(|| self.children.get(key).cloned())
    }

    pub fn generate_combinations(&self, family: &Family) -> Vec<Iht> {
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

    pub fn legend(&self) -> String {
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
    pub fn collapse_to_string(&self) -> String {
        // Helper function that chooses the separator:
        // if either allele is '.', use "|" regardless of the provided default separator.
        fn format_alleles(hap_a: &char, hap_b: &char, default_sep: &str) -> String {
            // Compare the alleles directly as chars.
            let sep = if *hap_a == '.' || *hap_b == '.' {
                "|"
            } else {
                default_sep
            };
            format!("{}{}{} ", hap_a, sep, hap_b)
        }

        // Sort and add founders to the string.
        let mut sorted_founders: Vec<_> = self.founders.iter().collect();
        sorted_founders.sort_by_key(|(id, _)| *id);
        let founders_str: String = sorted_founders
            .iter()
            .map(|(_, (hap_a, hap_b))| format_alleles(hap_a, hap_b, "/"))
            .collect::<Vec<_>>()
            .join(" ");

        // Sort and add children to the string.
        let mut sorted_children: Vec<_> = self.children.iter().collect();
        sorted_children.sort_by_key(|(id, _)| *id);
        let children_str: String = sorted_children
            .iter()
            .map(|(_, (hap_a, hap_b))| format_alleles(hap_a, hap_b, "|"))
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
    pub fn merge(&mut self, other: &Iht) {
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
    /// Assigns genotype alleles to both founders and children based on founder allele matches.
    pub fn assign_genotypes(
        &self,
        founder_genotypes: &HashMap<String, Vec<GenotypeAllele>>,
        sort_alleles: bool,
    ) -> (
        HashMap<char, GenotypeAllele>,
        HashMap<String, Vec<GenotypeAllele>>,
    ) {
        let mut updated_genotypes = HashMap::new();
        let mut founder_allele_map: HashMap<char, GenotypeAllele> = HashMap::new();

        // Step 1: Build `founder_allele_map` from founder genotypes
        for (founder, (allele1, allele2)) in &self.founders {
            if let Some(founder_alleles) = founder_genotypes.get(founder) {
                /*
                println!(
                    "fu {:?} {:?} a1:{} a2:{}",
                    founder, founder_alleles, allele1, allele2
                );
                */
                if founder_alleles.len() >= 2 {
                    founder_allele_map.insert(*allele1, founder_alleles[0]);
                    founder_allele_map.insert(*allele2, founder_alleles[1]);
                }
            }
        }

        //println!(" {:?}", founder_allele_map);

        // Step 2: Assign genotypes based on `founder_allele_map`
        for (individual, (hap1, hap2)) in self.founders.iter().chain(self.children.iter()) {
            let allele1 = founder_allele_map
                .get(hap1)
                .copied()
                .unwrap_or(GenotypeAllele::UnphasedMissing);
            let allele2 = founder_allele_map
                .get(hap2)
                .copied()
                .unwrap_or(GenotypeAllele::UnphasedMissing);

            // Ensure consistent ordering
            let sorted_genotypes = if allele1.index() > allele2.index() && sort_alleles {
                vec![allele2, allele1]
            } else {
                vec![allele1, allele2]
            };

            updated_genotypes.insert(individual.clone(), sorted_genotypes);
        }

        (founder_allele_map, updated_genotypes)
    }
    /// Generate all possible founder allele orientations
    pub fn founder_phase_orientations(&self) -> Vec<Iht> {
        let mut results = Vec::new();

        // Collect founder IDs for iteration
        let founder_ids: Vec<String> = self.founders.keys().cloned().collect();

        // Generate all combinations of allele swaps
        let num_founders = founder_ids.len();
        let all_swaps = (0..(1 << num_founders)) // 2^num_founders combinations
            .map(|i| {
                founder_ids
                    .iter()
                    .enumerate()
                    .map(|(j, id)| {
                        let (allele_a, allele_b) = self.founders[id];
                        if (i >> j) & 1 == 1 {
                            (id.clone(), (allele_b, allele_a)) // Swap alleles
                        } else {
                            (id.clone(), (allele_a, allele_b)) // Keep original
                        }
                    })
                    .collect::<HashMap<String, (char, char)>>()
            });

        // Create new `Iht` instances for each orientation
        for swapped_founders in all_swaps {
            let mut new_iht = self.clone();
            new_iht.founders = swapped_founders;
            results.push(new_iht);
        }

        results
    }

    /// Returns a HashSet of all non-'?' characters present in the children.
    pub fn get_non_missing_child_alleles(&self) -> HashSet<char> {
        let mut non_missing: HashSet<char> = HashSet::new();

        for (_, (hap_a, hap_b)) in &self.children {
            if *hap_a != '?' {
                non_missing.insert(*hap_a);
            }
            if *hap_b != '?' {
                non_missing.insert(*hap_b);
            }
        }

        non_missing
    }

    /// Returns a HashSet of all non-'?' characters present in the children
    /// **only for children of founders with multiple children**.
    pub fn get_flipable_alleles(&self, family: &Family) -> HashSet<char> {
        let mut non_missing: HashSet<char> = HashSet::new();

        // Identify founders with multiple children
        let founders_with_multiple_children: HashSet<String> = family
            .get_founders_with_multiple_children()
            .into_iter()
            .map(|f| f.id())
            .collect();

        // Iterate through children and only collect alleles for those with multi-child founders
        for (child_id, (hap_a, hap_b)) in &self.children {
            if let Some((father, mother)) = family.get_parents(child_id) {
                if founders_with_multiple_children.contains(&father)
                    || founders_with_multiple_children.contains(&mother)
                {
                    if *hap_a != '?' {
                        non_missing.insert(*hap_a);
                    }
                    if *hap_b != '?' {
                        non_missing.insert(*hap_b);
                    }
                }
            }
        }

        non_missing
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

/// Parses the IHT file and loads it into a Vec<IhtVec>
/// Supports both "|" and "/" as valid field separators.
pub fn parse_ihtv2_file(file_path: &str, founder_count: usize) -> Vec<IhtVec> {
    let file = File::open(file_path).expect("Failed to open the file");
    let reader = BufReader::new(file);

    let mut iht_vecs = Vec::new();
    let mut lines = reader.lines();

    // Read header
    let header_line = lines
        .next()
        .expect("Missing header")
        .expect("Invalid UTF-8");
    let headers: Vec<&str> = header_line.split_whitespace().collect();

    // Extract column indices for individuals
    let individual_start_idx = 3; // The first three columns are chrom, start, end
    let marker_count_idx = headers.len() - 2;
    let length_idx = headers.len() - 1;

    // Read data lines
    for line in lines {
        let line = line.expect("Invalid UTF-8");
        let fields: Vec<&str> = line.split_whitespace().collect();

        if fields.len() < headers.len() {
            continue; // Skip malformed lines
        }

        let chrom = fields[0].to_string();
        let start: i64 = fields[1].parse().expect("Invalid start position");
        let end: i64 = fields[2].parse().expect("Invalid end position");

        let marker_count: usize = fields[marker_count_idx]
            .parse()
            .expect("Invalid marker count");
        let length: usize = fields[length_idx].parse().expect("Invalid length");

        // Parse founders and children
        let mut founders = HashMap::new();
        let mut children = HashMap::new();

        let child_names = &headers[individual_start_idx + founder_count..marker_count_idx];

        for (i, &alleles) in fields[individual_start_idx..marker_count_idx]
            .iter()
            .enumerate()
        {
            // Split the field using both "|" and "/"
            let parts: Vec<&str> = alleles.split(|c| c == '|' || c == '/').collect();

            if parts.len() != 2 || parts[0].is_empty() || parts[1].is_empty() {
                continue; // Skip malformed entries
            }

            let individual_name = headers[individual_start_idx + i].to_string();
            let entry = (
                parts[0].chars().next().unwrap_or('?'), // Use '?' if missing
                parts[1].chars().next().unwrap_or('?'),
            );

            if i < founder_count {
                founders.insert(individual_name, entry);
            } else {
                children.insert(individual_name, entry);
            }
        }

        // Store in IhtVec struct
        iht_vecs.push(IhtVec {
            bed: BedRecord { chrom, start, end },
            iht: Iht { founders, children },
            non_missing_counts: HashMap::new(),
            count: marker_count,
        });
    }

    iht_vecs
}
