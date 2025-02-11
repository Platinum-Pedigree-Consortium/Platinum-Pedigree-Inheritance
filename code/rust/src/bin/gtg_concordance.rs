use clap::Parser;
use concordance::bed::BedRecord;
use concordance::iht::Iht;
use concordance::iht::IhtVec;
use concordance::utils::has_missing_alleles;

use concordance::iht::parse_ihtv2_file;
use concordance::ped::Family;
use concordance::utils::get_sample_depths;
use concordance::utils::is_vcf_indexed;
use log::debug;
use log::info;
use log::warn;
use log::LevelFilter;
use std::collections::HashMap;
use std::io;
use std::process;
use std::str;

use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{IndexedReader, Read};

/// Filter and Phase regions in haplotype map.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = " A tool filter/phase based on a haplotype (inheritance vectors) map."
)]
struct Args {
    /// A PED file documenting familial relationships
    #[arg(short, long)]
    ped: String,

    /// VCF file
    #[arg(long)]
    vcf: String,

    /// Inheritance file
    #[arg(short, long)]
    inheritance: String,

    /// Output prefix
    #[arg(long)]
    prefix: String,

    /// Minimum variant quality score
    #[arg(short, long, default_value_t = 20.0)]
    qual: f32,

    /// Minimum depth for all family members
    #[arg(short, long, default_value_t = 10)]
    depth: i32,

    /// Verbosity
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
}

/// Converts a `HashMap<String, (i32, Vec<GenotypeAllele>)>`
/// into a `HashMap<String, Vec<GenotypeAllele>>` by stripping out the depth.
fn convert_genotype_map(
    input: &HashMap<String, (i32, Vec<GenotypeAllele>)>,
) -> HashMap<String, Vec<GenotypeAllele>> {
    input
        .iter()
        .map(|(key, (_, alleles))| (key.clone(), alleles.clone())) // Remove depth, keep only alleles
        .collect()
}

/// Reads all records in a specific chromosome from a VCF file and extracts genotype alleles.
/// Converts the chromosome name to a sequence index before fetching records.
/// Returns a vector where each element is a tuple of:
/// - A `BedRecord` containing the position.
/// - A HashMap where keys are sample IDs and values are `(depth, alleles)`.
fn parse_vcf(
    reader: &mut IndexedReader,
    region: &BedRecord,
    minqual: f32,
) -> io::Result<Vec<(BedRecord, HashMap<String, (i32, Vec<GenotypeAllele>)>)>> {
    // Convert chromosome name to sequence index
    let chrom_id = match reader.header().name2rid(region.chrom.as_bytes()) {
        Ok(id) => id,
        Err(_) => {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "Chromosome not found",
            ))
        }
    };
    // Convert start and end positions safely
    let start: u64 = region.start.try_into().unwrap_or(0); // Ensure start is valid
    let end: Option<u64> = if region.end >= 0 {
        Some(region.end.try_into().unwrap_or(u64::MAX)) // Ensure end is valid
    } else {
        None // No end limit
    };

    // Fetch the records using the sequence index and converted start/end
    let rv = reader.fetch(chrom_id, start, end);

    // Initialize a vector to store genotype data for each site
    let mut records_genotype_map: Vec<(BedRecord, HashMap<String, (i32, Vec<GenotypeAllele>)>)> =
        Vec::new();

    if rv.is_err() {
        return Ok(records_genotype_map); // Return empty if fetch fails
    }

    // Extract sample names from the VCF header
    let header = reader.header();
    let samples: Vec<String> = header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect();

    // Iterate through all records in the specified chromosome region
    for result in reader.records() {
        let record = match result {
            Ok(rec) => rec,
            Err(_) => continue, // Skip invalid records
        };

        if record.qual() < minqual {
            continue;
        }

        let depths = match get_sample_depths(&record, &samples) {
            Some(d) => d,
            None => continue,
        };

        let genotypes = match record.genotypes() {
            Ok(g) => g,
            Err(_) => continue,
        };

        let mut record_map: HashMap<String, (i32, Vec<GenotypeAllele>)> = HashMap::new();

        // Extract genotype alleles for each sample
        for (i, sample) in samples.iter().enumerate() {
            let geno: rust_htslib::bcf::record::Genotype = genotypes.get(i);
            let mut alleles: Vec<GenotypeAllele> = geno.iter().cloned().collect();
            alleles.sort_by_key(|a| a.index().unwrap_or(i32::MAX.try_into().unwrap()));

            record_map.insert(sample.clone(), (*depths.get(sample).unwrap_or(&0), alleles));
        }

        // Store results
        records_genotype_map.push((
            BedRecord {
                chrom: region.chrom.clone(),
                start: record.pos() as i64,
                end: record.pos() as i64,
            },
            record_map,
        ));
    }

    Ok(records_genotype_map)
}

/// Prints two `HashMap<String, Vec<GenotypeAllele>>` side by side for comparison,
/// adding columns for mismatches (using `genotypes_match` logic) and Iht characters.
fn print_genotype_maps(
    map1: &HashMap<String, Vec<GenotypeAllele>>,
    map2: &HashMap<String, Vec<GenotypeAllele>>,
    iht: &Iht, // Iht structure containing founder and child allele characters
) {
    // Get a sorted list of all unique keys from both maps
    let mut all_keys: Vec<String> = map1.keys().chain(map2.keys()).cloned().collect();
    all_keys.sort();
    all_keys.dedup(); // Remove duplicates

    // Compute mismatch list using `genotypes_match`
    let mut mismatches = Vec::new();
    for (id, alleles1) in map1 {
        if let Some(alleles2) = map2.get(id) {
            if !genotypes_match(alleles1, alleles2) {
                mismatches.push(id.clone());
            }
        }
    }

    // Print the header
    println!(
        "{:<15} {:<25} {:<25} {:<10} {:<10}",
        "Sample", "Seen", "Expected", "Matching", "Iht"
    );
    println!("{}", "-".repeat(100));

    // Loop through all keys and print side-by-side values
    for key in all_keys {
        let alleles1 = map1
            .get(&key)
            .map_or_else(|| "None".to_string(), |v| format!("{:?}", v));
        let alleles2 = map2
            .get(&key)
            .map_or_else(|| "None".to_string(), |v| format!("{:?}", v));

        // Use precomputed mismatch list
        let mismatch_marker = if mismatches.contains(&key) {
            "YES"
        } else {
            "NO"
        };

        // Retrieve Iht characters
        let iht_chars = iht
            .get_alleles(&key)
            .map_or_else(|| "--".to_string(), |(a, b)| format!("{}{}", a, b));

        println!(
            "{:<15} {:<25} {:<25} {:<10} {:<10}",
            key, alleles1, alleles2, mismatch_marker, iht_chars
        );
    }
}

/// Compare two genotype maps and return a vector of IDs that don't match.
/// Ignores both `GenotypeAllele::UnphasedMissing` and `GenotypeAllele::PhasedMissing`.
fn compare_genotype_maps(
    map1: &HashMap<String, Vec<GenotypeAllele>>,
    map2: &HashMap<String, Vec<GenotypeAllele>>,
) -> Vec<String> {
    let mut mismatches = Vec::new();

    for (id, alleles1) in map1 {
        if let Some(alleles2) = map2.get(id) {
            if !genotypes_match(alleles1, alleles2) {
                mismatches.push(id.clone());
            }
        }
    }

    mismatches
}

/// Helper function to compare two genotype vectors while ignoring missing values
fn genotypes_match(alleles1: &Vec<GenotypeAllele>, alleles2: &Vec<GenotypeAllele>) -> bool {
    let filtered1: Vec<_> = alleles1.iter().filter(|a| !is_ignored_missing(a)).collect();
    let filtered2: Vec<_> = alleles2.iter().filter(|a| !is_ignored_missing(a)).collect();

    // If both filtered vectors are empty, consider them as matching
    if filtered1.is_empty() && filtered2.is_empty() {
        return true;
    }

    // Compare filtered vectors
    filtered1 == filtered2
}

/// Helper function to ignore both `UnphasedMissing` and `PhasedMissing`
fn is_ignored_missing(allele: &GenotypeAllele) -> bool {
    matches!(
        allele,
        GenotypeAllele::UnphasedMissing | GenotypeAllele::PhasedMissing
    )
}

fn find_best_phase_orientation(
    genotypes: &HashMap<String, (i32, Vec<GenotypeAllele>)>,
    iht_vec: &IhtVec,
) -> (Option<Iht>, Vec<String>) {
    let founder_phases = iht_vec.iht.founder_phase_orientations();
    let mut best_phase: Option<Iht> = None;
    let mut lowest_mismatch = usize::MAX;
    let mut best_mismatch_vec = Vec::new();

    // Convert genotypes once (no outer loop)
    let converted_genotypes = convert_genotype_map(genotypes);

    for phase in &founder_phases {
        println!("\n{}", "+".repeat(100));
        println!("{}", phase);
        let new_genotypes = phase.assign_genotypes(&converted_genotypes);
        let mismatch_vec = compare_genotype_maps(&converted_genotypes, &new_genotypes.1);
        let mismatch_count = mismatch_vec.len();

        print_genotype_maps(&converted_genotypes, &new_genotypes.1, phase);

        println!(
            "\nTotal Mismatch Count:{} Samples:{}\n",
            mismatch_count,
            mismatch_vec.join(","),
        );

        let mut sorted_keys: Vec<char> = new_genotypes.0.keys().cloned().collect();
        sorted_keys.sort(); // Sorting characters naturally

        // Print each key-value pair in sorted order
        for key in sorted_keys {
            println!(
                "iht2allele {} -> {:?}",
                key,
                new_genotypes.0.get(&key).unwrap()
            );
        }

        // Store the phase with the lowest mismatch
        if mismatch_count < lowest_mismatch {
            lowest_mismatch = mismatch_count;
            best_phase = Some(phase.clone());
            best_mismatch_vec = mismatch_vec.clone();
        }
    }
    println!("\n{}", "+".repeat(100));
    (best_phase, best_mismatch_vec)
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

    let family = Family::parse_ped_file(&args.ped).unwrap();

    let mut reader: IndexedReader =
        IndexedReader::from_path(&args.vcf).expect("Failure to read VCF file.");

    let mut output_vcf = args.prefix.clone();
    output_vcf += ".vcf";

    let iht_info = parse_ihtv2_file(&args.inheritance, family.founders().len());

    for v in iht_info {
        debug!("{} {} {} {}", v.bed.chrom, v.bed.start, v.bed.end, v.iht,);

        let genotypes = parse_vcf(&mut reader, &v.bed, args.qual).unwrap();
        let founder_phases = v.iht.founder_phase_orientations();

        let mut passing_count = 0;
        let mut failing_count = 0;
        let mut nocall_count = 0;

        for records in &genotypes {
            if has_missing_alleles(&records.1, 5) {
                nocall_count += 1;
                continue;
            }

            let best_results = find_best_phase_orientation(&records.1, &v);
            if !best_results.1.is_empty() {
                let mut issues = best_results.1;
                issues.sort();
                let converted_genos = convert_genotype_map(&records.1);
                let expected = best_results
                    .0
                    .as_ref()
                    .unwrap()
                    .assign_genotypes(&converted_genos);
                println!(
                    "\n\nBR {} {} {}",
                    records.0.start,
                    issues.len(),
                    issues.join(","),
                );
                print_genotype_maps(&converted_genos, &expected.1, &best_results.0.unwrap());

                failing_count += 1;
            } else {
                passing_count += 1;
            }
        }

        info!(
            "region {}:{}-{} has {} variants",
            v.bed.chrom,
            v.bed.start,
            v.bed.end,
            genotypes.len()
        );
        info!(
            "Passing {} Failing {} Nocall {}",
            passing_count, failing_count, nocall_count
        );
    }
}
