use log::{debug, warn};
use std::collections::HashMap;
use std::collections::VecDeque;
use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::path::Path;

#[derive(Debug)]
pub struct Individual {
    family_id: String,
    id: String,
    father_id: Option<String>,
    mother_id: Option<String>,
    sex: Option<u8>,
    phenotype: Option<u8>,
    children: Vec<String>,
}

impl Individual {
    pub fn new(
        family_id: String,
        id: String,
        father_id: Option<String>,
        mother_id: Option<String>,
        sex: Option<u8>,
        phenotype: Option<u8>,
    ) -> Self {
        Self {
            family_id,
            id,
            father_id,
            mother_id,
            sex,
            phenotype,
            children: Vec::new(),
        }
    }

    pub fn id(&self) -> String {
        self.id.clone()
    }

    pub fn no_parents(&self) -> bool {
        self.father_id.is_none() && self.mother_id.is_none()
    }

    pub fn has_children(&self) -> bool {
        !self.children.is_empty()
    }
    pub fn get_sex(&self) -> Option<u8> {
        self.sex
    }
    pub fn get_father_id(&self) -> Option<String> {
        self.father_id.clone()
    }

    pub fn get_mother_id(&self) -> Option<String> {
        self.mother_id.clone()
    }
}

pub struct Family {
    individuals: HashMap<String, Individual>,
}

impl Family {
    fn new(mut individuals: HashMap<String, Individual>) -> Self {
        let ids: Vec<String> = individuals.keys().cloned().collect();
        for id in ids {
            // Retrieve the individual and destructure its parent IDs
            let (father_id, mother_id) = if let Some(ind) = individuals.get(&id) {
                (ind.father_id.clone(), ind.mother_id.clone())
            } else {
                continue;
            };

            // Add the child ID to the father's children list
            if let Some(father_id) = father_id {
                if let Some(father) = individuals.get_mut(&father_id) {
                    father.children.push(id.clone());
                }
            }

            // Add the child ID to the mother's children list
            if let Some(mother_id) = mother_id {
                if let Some(mother) = individuals.get_mut(&mother_id) {
                    mother.children.push(id.clone());
                }
            }
        }
        Family { individuals }
    }

    pub fn len(&self) -> usize {
        self.individuals.len()
    }

    pub fn founders(&self) -> Vec<&Individual> {
        let mut founders: Vec<&Individual> = self
            .individuals
            .values()
            .filter(|ind| ind.no_parents() && ind.has_children())
            .collect();

        // Sort by ID (assuming Individual has a method `id()` that returns &str)
        founders.sort_by_key(|ind| ind.id());

        founders
    }

    /// Returns a list of offspring (individuals with parents) sorted by increasing depth.
    ///
    /// # Description
    /// - This function retrieves all individuals who are not founders (i.e., they have at least one parent).
    /// - It sorts the offspring based on their depth in the pedigree, ensuring that children appear before their descendants.
    ///
    /// # Returns
    /// - A `Vec<&Individual>` containing offspring sorted by increasing depth.
    ///
    /// # Example
    /// ```rust
    /// let offspring = family.offspring();
    /// for ind in offspring {
    ///     println!("ID: {}, Depth: {:?}", ind.id, family.get_individual_depths().iter().find(|(id, _)| id == &ind.id).unwrap().1);
    /// }
    /// ```
    pub fn offspring(&self) -> Vec<&Individual> {
        let depths = self.get_individual_depths(); // Get individuals sorted by depth

        depths
            .iter()
            .filter_map(|(id, _)| self.individuals.get(id)) // Get Individual references
            .filter(|ind| !ind.no_parents()) // Keep only non-founders
            .collect()
    }
    /// Returns indviduals who are parents of child
    pub fn get_parents(&self, child_id: &String) -> Option<(String, String)> {
        self.individuals.get(child_id).and_then(|child| {
            match (&child.father_id, &child.mother_id) {
                (Some(father), Some(mother)) => Some((father.clone(), mother.clone())),
                _ => None,
            }
        })
    }
    /// Returns the children of a given parent ID.
    pub fn get_children(&self, parent_id: &String) -> Vec<&String> {
        self.individuals
            .values()
            .filter(|ind| {
                ind.father_id.as_ref() == Some(parent_id)
                    || ind.mother_id.as_ref() == Some(parent_id)
            })
            .map(|ind| &ind.id)
            .collect()
    }
    /// Returns a vector of tuples `(sample_id, depth)` for every individual in the family.
    /// Depth is the distance from a founder.
    pub fn get_individual_depths(&self) -> Vec<(String, usize)> {
        let mut depths = HashMap::new(); // Track the maximum depth of each individual
        let mut queue: VecDeque<(String, usize)> = VecDeque::new();

        // Find all founders (individuals with no parents)
        let founders: Vec<&String> = self
            .individuals
            .iter()
            .filter_map(|(id, ind)| if ind.no_parents() { Some(id) } else { None })
            .collect();

        // Initialize the queue with founders
        for founder in founders {
            queue.push_back((founder.clone(), 0)); // Founders have depth 0
        }

        // Perform a breadth-first traversal
        while let Some((current_id, current_depth)) = queue.pop_front() {
            // Update the depth of the current individual to the maximum depth seen so far
            let entry = depths.entry(current_id.clone()).or_insert(current_depth);
            if current_depth > *entry {
                *entry = current_depth;
            }

            // Add children of the current individual to the queue with depth + 1
            if let Some(ind) = self.individuals.get(&current_id) {
                for child_id in &ind.children {
                    queue.push_back((child_id.clone(), current_depth + 1));
                }
            }
        }

        // Convert the HashMap to a sorted Vec<(String, usize)>
        let mut sorted_depths: Vec<(String, usize)> = depths.into_iter().collect();
        sorted_depths.sort_by(|a, b| a.1.cmp(&b.1).then(a.0.cmp(&b.0)));

        sorted_depths
    }
    /// Parses a PED file and constructs a `Family` instance.
    ///
    /// # Parameters:
    /// - `file_path`: Path to the PED file.
    ///
    /// # Returns:
    /// - A `Result` containing a `Family` object or an `io::Error` if parsing fails.
    pub fn parse_ped_file(file_path: &str) -> io::Result<Self> {
        let path = Path::new(file_path);
        let file = File::open(&path)?;
        let reader = BufReader::new(file);

        let mut individuals = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            let fields: Vec<&str> = line.split_whitespace().collect();

            if fields.len() < 6 {
                eprintln!("Skipping invalid line: {}", line);
                continue;
            }

            let family_id = fields[0].to_string();
            let id = fields[1].to_string();
            let father_id = if fields[2] == "NA" {
                None
            } else {
                Some(fields[2].to_string())
            };
            let mother_id = if fields[3] == "NA" {
                None
            } else {
                Some(fields[3].to_string())
            };
            let sex = fields[4].parse::<u8>().ok();
            let phenotype = fields[5].parse::<u8>().ok();

            let individual =
                Individual::new(family_id, id.clone(), father_id, mother_id, sex, phenotype);

            individuals.insert(id, individual);
        }

        if individuals.is_empty() {
            warn!("PED file was empty or malformed, no data loaded.");
        }
        debug!("Parsed pedigree okay");

        Ok(Family::new(individuals))
    }

    /// Returns a sorted vector of all individual IDs in the family.
    pub fn get_individuals_ids(&self) -> Vec<String> {
        let mut ids: Vec<String> = self.individuals.keys().cloned().collect();
        ids.sort(); // Sort the IDs in ascending order
        ids
    }

    /// Retrieves an `Individual` from the `Family` by their ID.
    ///
    /// # Arguments
    /// - `individual_id`: The ID of the individual to retrieve.
    ///
    /// # Returns
    /// - `Some(&Individual)` if the individual exists.
    /// - `None` if the individual is not found.
    pub fn get_individual(&self, individual_id: &str) -> Option<&Individual> {
        self.individuals.get(individual_id)
    }
}
