use clap::Parser;
use log::info;
use std::collections::HashSet;
use std::fmt;
use std::fs::read_to_string;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;

use log::LevelFilter;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    fasta: String,
    #[arg(short, long)]
    prefix: String,
    #[arg(short = 'v', long = "verbose", action = clap::ArgAction::Count, default_value_t = 0)]
    verbosity: u8,
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

    let file = File::open(args.fasta).unwrap();
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let lc = line.as_ref().unwrap().as_bytes()[0];

        if lc == b'>' {
            let lc_parts = line.unwrap().replace(">", "");
            println!(">{}_{}", args.prefix, lc_parts);
        } else {
            println!("{}", line.unwrap());
        }
    }
}
