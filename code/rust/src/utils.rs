use log::warn;
use rust_htslib::bcf::record::GenotypeAllele;
use rust_htslib::bcf::IndexedReader;
use rust_htslib::bcf::Read;
use rust_htslib::bcf::Record;
use std::collections::HashMap;
use std::io;
use std::path::Path;

pub fn get_sample_depths(record: &Record, samples: &Vec<String>) -> Option<HashMap<String, i32>> {
    // Try DP field first
    if let Ok(dp_values) = record.format(b"DP").integer() {
        let mut depths = HashMap::new();
        for (i, dp_sample) in dp_values.iter().enumerate() {
            if let Some(&depth) = dp_sample.get(0) {
                if let Some(sample_name) = samples.get(i) {
                    depths.insert(sample_name.clone(), depth);
                }
            }
        }
        return Some(depths);
    }

    // Fallback to AD field sum
    if let Ok(ad_values) = record.format(b"AD").integer() {
        let mut depths = HashMap::new();
        for (i, ad_sample) in ad_values.iter().enumerate() {
            let sum: i32 = ad_sample.iter().copied().sum();
            if let Some(sample_name) = samples.get(i) {
                depths.insert(sample_name.clone(), sum);
            }
        }
        return Some(depths);
    }

    // Fallback to SD field sum
    if let Ok(sd_values) = record.format(b"SD").integer() {
        let mut depths = HashMap::new();
        for (i, sd_sample) in sd_values.iter().enumerate() {
            let sum: i32 = sd_sample.iter().copied().sum();
            if let Some(sample_name) = samples.get(i) {
                depths.insert(sample_name.clone(), sum);
            }
        }
        return Some(depths);
    }

    // Neither DP, AD, nor SD present
    warn!("No field DP, AD, or SD, skipping variant.");
    None
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

pub fn unphase_genotype(alleles: &[GenotypeAllele]) -> Vec<GenotypeAllele> {
    let mut sorted_alleles: Vec<GenotypeAllele> = alleles
        .iter()
        .map(|a| match a {
            GenotypeAllele::PhasedMissing => GenotypeAllele::UnphasedMissing, // Convert phased missing to unphased missing
            GenotypeAllele::Phased(x) => GenotypeAllele::Unphased(*x), // Convert phased alleles to unphased
            _ => *a,                                                   // Keep everything else as is
        })
        .collect();

    // Sort by allele index, ensuring missing values come first
    sorted_alleles.sort_by_key(|a| match a {
        GenotypeAllele::UnphasedMissing => u32::MIN, // Ensure missing calls come first
        GenotypeAllele::Unphased(x) => *x as u32,    // Sort alleles numerically
        _ => u32::MAX,                               // Shouldn't happen, but keeps safety
    });

    sorted_alleles
}
#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bcf::record::GenotypeAllele;

    #[test]
    fn test_unphase_homozygous() {
        let alleles = vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(1)]
        );
    }

    #[test]
    fn test_unphase_heterozygous() {
        let alleles = vec![GenotypeAllele::Unphased(1), GenotypeAllele::Unphased(0)];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]
        );
    }

    #[test]
    fn test_unphase_phased_heterozygous() {
        let alleles = vec![GenotypeAllele::Phased(1), GenotypeAllele::Phased(0)];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(1)]
        );
    }

    #[test]
    fn test_unphase_phased_homozygous() {
        let alleles = vec![GenotypeAllele::Phased(2), GenotypeAllele::Phased(2)];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![GenotypeAllele::Unphased(2), GenotypeAllele::Unphased(2)]
        );
    }

    #[test]
    fn test_unphase_with_nocall() {
        let alleles = vec![GenotypeAllele::Unphased(1), GenotypeAllele::UnphasedMissing];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![GenotypeAllele::UnphasedMissing, GenotypeAllele::Unphased(1)]
        );
    }

    #[test]
    fn test_unphase_phased_with_nocall() {
        let alleles = vec![GenotypeAllele::Phased(1), GenotypeAllele::PhasedMissing];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![GenotypeAllele::UnphasedMissing, GenotypeAllele::Unphased(1)]
        );
    }

    #[test]
    fn test_unphase_nocall_only() {
        let alleles = vec![
            GenotypeAllele::UnphasedMissing,
            GenotypeAllele::UnphasedMissing,
        ];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::UnphasedMissing
            ]
        );
    }

    #[test]
    fn test_unphase_multi_allelic() {
        let alleles = vec![GenotypeAllele::Unphased(2), GenotypeAllele::Unphased(0)];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![GenotypeAllele::Unphased(0), GenotypeAllele::Unphased(2)]
        );
    }

    #[test]
    fn test_unphase_larger_genotype() {
        let alleles = vec![
            GenotypeAllele::Unphased(3),
            GenotypeAllele::Unphased(1),
            GenotypeAllele::Unphased(2),
            GenotypeAllele::Unphased(0),
        ];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![
                GenotypeAllele::Unphased(0),
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2),
                GenotypeAllele::Unphased(3)
            ]
        );
    }

    #[test]
    fn test_unphase_with_mixed_nocalls() {
        let alleles = vec![
            GenotypeAllele::Unphased(1),
            GenotypeAllele::UnphasedMissing,
            GenotypeAllele::Unphased(2),
        ];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2)
            ]
        );
    }

    #[test]
    fn test_unphase_phased_mixed_nocalls() {
        let alleles = vec![
            GenotypeAllele::Phased(1),
            GenotypeAllele::PhasedMissing,
            GenotypeAllele::Phased(2),
        ];
        let result = unphase_genotype(&alleles);
        assert_eq!(
            result,
            vec![
                GenotypeAllele::UnphasedMissing,
                GenotypeAllele::Unphased(1),
                GenotypeAllele::Unphased(2)
            ]
        );
    }
}
