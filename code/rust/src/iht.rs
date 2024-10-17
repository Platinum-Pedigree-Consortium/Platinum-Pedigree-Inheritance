use csv::ReaderBuilder;
use std::collections::HashMap;
use std::collections::HashSet;

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
        println!("{}", ihtblock);
        inht_info.push(ihtblock);
    }
    return inht_info;
}
