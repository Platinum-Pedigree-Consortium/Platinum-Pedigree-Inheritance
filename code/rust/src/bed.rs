use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

// Define a structure to represent a BED record
#[derive(Debug)]
pub struct BedRecord {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl BedRecord {
    // Function to parse a line from the BED file into a BedRecord
    pub fn from_line(line: &str) -> Result<BedRecord, &'static str> {
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 3 {
            return Err("Invalid BED format");
        }

        let chrom = fields[0].to_string();
        let start = fields[1]
            .parse::<u64>()
            .map_err(|_| "Invalid start position")?;
        let end = fields[2]
            .parse::<u64>()
            .map_err(|_| "Invalid end position")?;

        Ok(BedRecord { chrom, start, end })
    }
}

// Function to read a BED file
pub fn read_bed_file<P: AsRef<Path>>(file_path: P) -> io::Result<Vec<BedRecord>> {
    let file = File::open(file_path)?;
    let reader = io::BufReader::new(file);
    let mut records = Vec::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue; // Skip comment lines and empty lines
        }

        match BedRecord::from_line(&line) {
            Ok(record) => records.push(record),
            Err(e) => eprintln!("Skipping line due to error: {}", e),
        }
    }

    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed_record_from_line() {
        let line = "chr1 1000 2000";
        let record = BedRecord::from_line(line).unwrap();
        assert_eq!(record.chrom, "chr1");
        assert_eq!(record.start, 1000);
        assert_eq!(record.end, 2000);
    }

    #[test]
    fn test_read_bed_file() {
        let bed_data = "chr1 1000 2000\nchr2 1500 2500\n";
        let temp_file_path = "test.bed";
        std::fs::write(temp_file_path, bed_data).unwrap();

        let records = read_bed_file(temp_file_path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].chrom, "chr1");
        assert_eq!(records[1].chrom, "chr2");

        std::fs::remove_file(temp_file_path).unwrap();
    }
}
