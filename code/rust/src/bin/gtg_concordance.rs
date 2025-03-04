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
use std::fs::OpenOptions;
use std::io;
use std::io::Write;
use std::path::Path;
use std::process;
use std::str;

use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::{Format, Header, IndexedReader, Read, Writer};

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
    #[arg(short, long, default_value_t = 5)]
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

/// Extracts genotype data from a single VCF record.
fn parse_vcf_record(
    chrom: &String,
    record: &rust_htslib::bcf::Record,
    samples: &Vec<String>,
) -> (BedRecord, HashMap<String, (i32, Vec<GenotypeAllele>)>) {
    let mut record_map: HashMap<String, (i32, Vec<GenotypeAllele>)> = HashMap::new();

    if let Some(depths) = get_sample_depths(record, samples) {
        if let Ok(genotypes) = record.genotypes() {
            for (i, sample) in samples.iter().enumerate() {
                let geno: rust_htslib::bcf::record::Genotype = genotypes.get(i);
                let mut alleles: Vec<GenotypeAllele> = geno.iter().cloned().collect();
                alleles.sort_by_key(|a| a.index().unwrap_or(i32::MAX.try_into().unwrap()));
                record_map.insert(sample.clone(), (*depths.get(sample).unwrap_or(&0), alleles));
            }
        }
    }

    (
        BedRecord {
            chrom: chrom.clone(),
            start: record.pos() as i64,
            end: record.pos() as i64,
        },
        record_map,
    )
}

/// Returns a formatted string comparing two `HashMap<String, Vec<GenotypeAllele>>`
/// side by side, adding columns for mismatches and Iht characters.
fn format_genotype_maps(
    map1: &HashMap<String, Vec<GenotypeAllele>>,
    map2: &HashMap<String, Vec<GenotypeAllele>>,
    iht: &Iht, // Iht structure containing founder and child allele characters
) -> String {
    let mut output = String::new();

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

    // Format the header
    output.push_str(&format!(
        "{:<15} {:<25} {:<25} {:<10} {:<10}\n",
        "Sample", "Seen", "Expected", "Mismatch", "Inheritance"
    ));
    output.push_str(&"-".repeat(100));
    output.push('\n');

    // Loop through all keys and format side-by-side values
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

        output.push_str(&format!(
            "{:<15} {:<25} {:<25} {:<10} {:<10}\n",
            key, alleles1, alleles2, mismatch_marker, iht_chars
        ));
    }

    output
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
    fam: &Family,
) -> (Option<Iht>, Vec<String>) {
    let founder_phases = iht_vec.iht.founder_phase_orientations();
    let mut best_phase: Option<Iht> = None;
    let mut lowest_mismatch = usize::MAX;
    let mut best_mismatch_vec = Vec::new();

    // Convert genotypes once (no outer loop)
    let converted_genotypes = convert_genotype_map(genotypes);

    for phase in &founder_phases {
        // println!("\n{}", "+".repeat(100));
        // println!("P: {}", phase.collapse_to_string());
        let new_genotypes = phase.assign_genotypes(&converted_genotypes, true);
        let mismatch_vec = compare_genotype_maps(&converted_genotypes, &new_genotypes.1);
        let mismatch_count = mismatch_vec.len();

        /*
        print_genotype_maps(&converted_genotypes, &new_genotypes.1, phase);

        println!(
            "\nTotal Mismatch Count:{} Samples:{}\n",
            mismatch_count,
            mismatch_vec.join(","),
        );
        */

        let mut sorted_keys: Vec<char> = new_genotypes.0.keys().cloned().collect();
        sorted_keys.sort(); // Sorting characters naturally

        // Print each key-value pair in sorted order
        /*
        for key in sorted_keys {
            println!(
                "iht2allele {} -> {:?}",
                key,
                new_genotypes.0.get(&key).unwrap()
            );
        }
        println!("\n");
        */

        // Store the phase with the lowest mismatch
        if mismatch_count < lowest_mismatch {
            lowest_mismatch = mismatch_count;
            best_phase = Some(phase.clone());
            best_mismatch_vec = mismatch_vec.clone();
        }
    }
    // println!("\n{}", "+".repeat(100));
    (best_phase, best_mismatch_vec)
}

fn calculate_fraction(numerator: f64, denominator: f64) -> Option<f64> {
    if denominator == 0.0 {
        None // Return None if division by zero
    } else {
        Some(numerator / denominator)
    }
}

fn write_vcfs(
    reader: &mut IndexedReader,
    passing: &mut Writer,
    failing: &mut Writer,
    region: &BedRecord,
) {
    // Convert chromosome name to sequence index
    let chrom_id = reader.header().name2rid(region.chrom.as_bytes()).unwrap();
    // Convert start and end positions safely
    let start: u64 = region.start.try_into().unwrap_or(0); // Ensure start is valid
    let end: Option<u64> = if region.end >= 0 {
        Some(region.end.try_into().unwrap_or(u64::MAX)) // Ensure end is valid
    } else {
        None // No end limit
    };

    // Fetch the records using the sequence index and converted start/end
    let rv = reader.fetch(chrom_id, start, end).unwrap();

    for result in reader.records() {
        let record = match result {
            Ok(rec) => rec,
            Err(_) => continue, // Skip invalid records
        };
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

    let header = reader.header().clone();
    let wheader: Header = Header::from_template(&header);

    // Extract sample names from the VCF header
    let header = reader.header();
    let samples: Vec<String> = header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect();

    let mut outpassingvcf = Writer::from_path(
        Path::new(&format!("{}.pass.vcf", args.prefix)),
        &wheader,
        true,
        Format::Vcf,
    )
    .unwrap();

    let mut outfailvcf = Writer::from_path(
        Path::new(&format!("{}.fail.vcf", args.prefix)),
        &wheader,
        true,
        Format::Vcf,
    )
    .unwrap();

    let mut output_stats_fn = args.prefix.clone();
    output_stats_fn += ".filtering_stats.txt";

    let mut failed_sites_fn = args.prefix.clone();
    failed_sites_fn += ".failed_sites.txt";

    // Open the file with write and create options
    let mut stats_fh = OpenOptions::new()
        .write(true) // Open for writing
        .create(true) // Create if it doesn't exist
        .open(output_stats_fn)
        .unwrap();

    // Open the file with write and create options
    let mut failed_fh = OpenOptions::new()
        .write(true) // Open for writing
        .create(true) // Create if it doesn't exist
        .open(failed_sites_fn)
        .unwrap();

    failed_fh
        .write(format!("chrom start sample n-samples iht-region\n").as_bytes())
        .unwrap();

    stats_fh
        .write(format!("#chrom start end passing failing nocall passing-rate one-off\n").as_bytes())
        .unwrap();

    let iht_info = parse_ihtv2_file(&args.inheritance, family.founders().len());

    let mut total_passing_count = 0;
    let mut total_failing_count = 0;
    let mut total_nocall_count = 0;
    let mut total_low_qual_count = 0;

    for v in iht_info {
        debug!("{} {} {} {}", v.bed.chrom, v.bed.start, v.bed.end, v.iht);

        let mut passing_count = 0;
        let mut failing_count = 0;
        let mut nocall_count = 0;
        let mut low_qual_count = 0;

        let mut failed_singletons: HashMap<String, i32> = HashMap::new();

        // Convert start and end positions safely
        let start: u64 = v.bed.start.try_into().unwrap_or(0); // Ensure start is valid
        let end: Option<u64> = if v.bed.end >= 0 {
            Some(v.bed.end.try_into().unwrap_or(u64::MAX)) // Ensure end is valid
        } else {
            None // No end limit
        };

        let chrom_id = reader.header().name2rid(v.bed.chrom.as_bytes()).unwrap();

        // Fetch the records using the sequence index and converted start/end
        let rv = reader.fetch(chrom_id, start, end);

        for r in reader.records() {
            let record = r.unwrap();
            let parsed_record = parse_vcf_record(&v.bed.chrom, &record, &samples);

            if record.qual() < args.qual {
                low_qual_count += 1;
                outfailvcf.write(&record);
                continue;
            }

            if has_missing_alleles(&parsed_record.1, args.depth) {
                nocall_count += 1;
                outfailvcf.write(&record);
                continue;
            }

            let best_results = find_best_phase_orientation(&parsed_record.1, &v, &family);
            if !best_results.1.is_empty() {
                let mut issues = best_results.1;
                issues.sort();

                for s in issues.iter() {
                    failed_fh
                        .write(
                            format!(
                                "{} {} {} {} {}-{}\n",
                                parsed_record.0.chrom,
                                parsed_record.0.start,
                                s,
                                issues.len(),
                                v.bed.start,
                                v.bed.end,
                            )
                            .as_bytes(),
                        )
                        .unwrap();
                }

                if issues.len() == 1 {
                    for i in &issues {
                        *failed_singletons.entry(i.clone()).or_insert(0) += 1;
                    }
                }

                let converted_genos = convert_genotype_map(&parsed_record.1);
                let expected = best_results
                    .0
                    .as_ref()
                    .unwrap()
                    .assign_genotypes(&converted_genos, true);

                if issues.len() == 1 {
                    let geno_table = format_genotype_maps(
                        &converted_genos,
                        &expected.1,
                        &best_results.0.unwrap(),
                    );
                    debug!(
                        "Genotype filtering info:\n
                        \n{}\n Position info -  {}:{}-{} {} {} {}\n",
                        geno_table,
                        v.bed.chrom,
                        v.bed.start,
                        v.bed.end,
                        parsed_record.0.start,
                        issues.len(),
                        issues.join(","),
                    );
                }
                outfailvcf.write(&record).unwrap();
                failing_count += 1;
            } else {
                let mut new_rec = record.clone();
                let mut new_gts: Vec<GenotypeAllele> = Vec::new();
                let convert_genos = convert_genotype_map(&parsed_record.1);
                let mut br = best_results.0.unwrap().clone();
                let loaded = br.assign_genotypes(&convert_genos, false);

                for s in &samples {
                    if loaded.1.contains_key(s) {
                        let alleles = loaded.1.get(s).unwrap();

                        new_gts.push(*alleles.get(0).unwrap());

                        if let Some(second_allele) = alleles.get(1).unwrap().index() {
                            new_gts.push(GenotypeAllele::Phased(second_allele as i32));
                        } else {
                            new_gts.push(*alleles.get(1).unwrap());
                        }
                    } else {
                        new_gts.push(GenotypeAllele::UnphasedMissing);
                        new_gts.push(GenotypeAllele::UnphasedMissing);
                    }
                }
                new_rec.push_genotypes(&new_gts).unwrap();
                outpassingvcf.write(&new_rec).unwrap();

                passing_count += 1;
            }
        }

        total_passing_count += passing_count;
        total_failing_count += failing_count;
        total_nocall_count += nocall_count;
        total_low_qual_count += low_qual_count;

        info!(
            "Region {}:{}-{} - Passing: {} Failing: {} Nocall: {} Low-qual: {}",
            v.bed.chrom,
            v.bed.start,
            v.bed.end,
            passing_count,
            failing_count,
            nocall_count,
            low_qual_count
        );

        let singletons = failed_singletons
            .iter()
            .map(|(key, value)| format!("{}:{}", key, value))
            .collect::<Vec<_>>()
            .join(";");

        stats_fh
            .write(
                format!(
                    "{} {} {} {} {} {} {:.3} {}\n",
                    v.bed.chrom,
                    v.bed.start,
                    v.bed.end,
                    passing_count,
                    failing_count,
                    nocall_count,
                    calculate_fraction(
                        passing_count as f64,
                        (passing_count + failing_count) as f64
                    )
                    .unwrap(),
                    singletons,
                )
                .as_bytes(),
            )
            .unwrap();
    }
    info!(
        "Total - Passing: {} Failing: {} Nocall {} Low-qual {}",
        total_passing_count, total_failing_count, total_nocall_count, total_low_qual_count
    );
}
