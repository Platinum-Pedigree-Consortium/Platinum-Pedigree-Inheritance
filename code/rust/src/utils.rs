use log::warn;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::IndexedReader;
use rust_htslib::bcf::Read;
use rust_htslib::bcf::Record;
use std::collections::HashMap;
use std::io;
use std::path::Path;

pub fn get_sample_depths(record: &Record, samples: &Vec<String>) -> Option<HashMap<String, i32>> {
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

pub fn is_vcf_indexed(vcf_path: &str) -> io::Result<bool> {
    let index_path = format!("{}.tbi", vcf_path);
    if Path::new(&index_path).exists() {
        Ok(true)
    } else {
        warn!("VCF file {} is not indexed.", vcf_path);
        Ok(false)
    }
}

/// Checks if any sample has missing alleles or a sequencing depth of zero.
///
/// # Arguments
/// - `genotype_data`: A `HashMap` mapping sample IDs to their depth and genotype alleles.
/// - `min depth` : A i32
/// # Returns
/// - `true` if any sample has a depth of less than 10 or missing.
/// - `false` if all samples have valid depths and alleles.
pub fn has_missing_alleles(
    genotype_data: &HashMap<String, (i32, Vec<GenotypeAllele>)>,
    min: i32,
) -> bool {
    genotype_data.values().any(|(depth, allele_vec)| {
        *depth < min || allele_vec.iter().any(|allele| allele.index().is_none())
    })
}

/// Computes the mean and standard deviation of a vector of values.
fn compute_mean_std(values: &[i32]) -> (f64, f64) {
    if values.is_empty() {
        return (0.0, 0.0);
    }

    let mean = values.iter().map(|&x| x as f64).sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|&x| (x as f64 - mean).powi(2))
        .sum::<f64>()
        / values.len() as f64;
    let std_dev = variance.sqrt();

    (mean, std_dev)
}

/// Reads a VCF file, extracts per-sample `DP`, and calculates mean/std for each sample.
/// Returns a HashMap where the key is the sample name and the value is a tuple (mean, std_dev).
pub fn extract_depth_statistics(
    reader: &mut IndexedReader,
    chrom: u32,
) -> io::Result<HashMap<String, (f64, f64)>> {
    // HashMap to store depth values per sample
    let mut sample_depths: HashMap<String, Vec<i32>> = HashMap::new();

    let mut sample_stats: HashMap<String, (f64, f64)> = HashMap::new();

    let rv = reader.fetch(chrom, 0, None);

    match rv {
        Err(_) => return Ok(sample_stats),
        _ => {}
    }

    // Retrieve sample names from the VCF header
    let header = reader.header();
    let sample_names: Vec<String> = header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect();

    for sample in &sample_names {
        sample_depths.insert(sample.clone(), Vec::new());
    }

    // Iterate through records in the VCF
    for record in reader.records() {
        let record = record.unwrap();

        // Extract per-sample `DP` values
        if let Ok(dp_values) = record.format(b"DP").integer() {
            for (i, sample) in sample_names.iter().enumerate() {
                if let Some(&dp) = dp_values.get(i).and_then(|dp_vec| dp_vec.first()) {
                    sample_depths.get_mut(sample).unwrap().push(dp);
                }
            }
        }
    }

    // Calculate statistics and store in a new HashMap

    for (sample, depths) in &sample_depths {
        let (mean, std_dev) = compute_mean_std(depths);
        sample_stats.insert(sample.clone(), (mean, std_dev));
    }

    Ok(sample_stats)
}
