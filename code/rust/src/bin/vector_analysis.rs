use clap::Parser;
use csv::ReaderBuilder;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::io::Read as IoRead;

/// quick script to calculate stats over the inheritance vectors
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// CSV containing the inheritance vectors
    #[arg(short, long)]
    inheritance: String,
}

fn main() {
    let args = Args::parse();
    let mut reader = ReaderBuilder::new().from_path(args.inheritance);

    let header = reader
        .as_mut()
        .expect("Error reading inheritance CSV header")
        .headers()
        .expect("Error reading inheritance CSV header")
        .clone();

    let mut countsum: HashMap<i32, i128> = HashMap::new();
    let mut hapsums: HashMap<char, i128> = HashMap::new();
    *hapsums.entry('A').or_insert(0);
    *hapsums.entry('B').or_insert(0);
    *hapsums.entry('C').or_insert(0);
    *hapsums.entry('D').or_insert(0);

    let mut hap_obs: HashMap<char, Vec<i32>> = HashMap::new();

    for record in reader.expect("Error reading inheritance CSV.").records() {
        let mut counter: HashMap<char, i32> = HashMap::new();
        *counter.entry('A').or_insert(0);
        *counter.entry('B').or_insert(0);
        *counter.entry('C').or_insert(0);
        *counter.entry('D').or_insert(0);

        let mut parts = record.as_ref().unwrap();

        let start: i128 = parts.get(1).unwrap().parse().unwrap();
        let end: i128 = parts.get(2).unwrap().parse().unwrap();
        let len = end - start;

        for i in 3..(parts.len() - 2) {
            let haps: Vec<char> = parts.get(i).unwrap().chars().collect();
            *counter.entry(haps[0]).or_insert(0) += 1;
            *counter.entry(haps[1]).or_insert(0) += 1;
            //  println!("p: {}", parts.get(i).unwrap());

            if !parts.get(0).unwrap().contains("chrX") {
                *hapsums.entry(haps[0]).or_insert(0) += len;
                *hapsums.entry(haps[1]).or_insert(0) += len;
            }
        }

        let mut out: String = "".to_string();
        for i in parts {
            out = format!("{} {}", out, i);
        }

        let mut ncov = 0;

        for i in &counter {
            if *i.1 > 0 {
                ncov += 1;
            }
        }

        *countsum.entry(ncov).or_insert(0) += len;

        let mut paint = "NA".to_string();
        if *counter.get(&'A').unwrap() == 0 {
            paint = "A".to_string();
        }
        if *counter.get(&'B').unwrap() == 0 {
            paint = "B".to_string();
        }
        if *counter.get(&'C').unwrap() == 0 {
            paint = "C".to_string();
        }
        if *counter.get(&'D').unwrap() == 0 {
            paint = "D".to_string();
        }

        out = format!(
            "{} {} {} {} {} {} {} {}",
            out,
            ncov,
            len,
            counter.get(&'A').unwrap(),
            counter.get(&'B').unwrap(),
            counter.get(&'C').unwrap(),
            counter.get(&'D').unwrap(),
            paint,
        );

        println!("{}", out);

        if !parts.get(0).unwrap().contains("chrX") {
            hap_obs
                .entry('A')
                .or_insert_with(Vec::new)
                .push(*counter.get(&'A').unwrap());
            hap_obs
                .entry('B')
                .or_insert_with(Vec::new)
                .push(*counter.get(&'B').unwrap());
            hap_obs
                .entry('C')
                .or_insert_with(Vec::new)
                .push(*counter.get(&'C').unwrap());
            hap_obs
                .entry('D')
                .or_insert_with(Vec::new)
                .push(*counter.get(&'D').unwrap());
        }
    }

    for l in &hapsums {
        println!("TO {} {}", l.0, l.1);
    }
    for i in hap_obs {
        for j in i.1 {
            println!("HS\t{}\t{}", i.0, j);
        }
    }
    for i in countsum {
        println!("CS\t{}\t{}", i.0, i.1);
    }
}
