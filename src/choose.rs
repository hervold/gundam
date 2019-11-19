extern crate gundam;
#[macro_use]
extern crate log;
extern crate bio;
extern crate env_logger;

use bio::pattern_matching::pssm::DNAMotif;
use gundam::*;
use std::env;
use std::process::exit;

fn main() {
    let _ = env_logger::init();
    let args = env::args().collect::<Vec<String>>();
    if args.len() != 3 {
        println!(
            "invalid # of args: {}.  usage: gundam <pos.fa> <neg.fa>",
            args.len()
        );
        exit(1);
    }

    let kmer_specs = vec![
        (KMER_LEN, MIN_GAP, MAX_GAP), // 5, 0, 20
        (4, 0, 2),
        (3, 0, 2),
    ];
    for (i, j, k, f) in
        DyadMotif::<DNAMotif>::passing_kmers(args[1].as_str(), args[2].as_str(), &kmer_specs)
    {
        println!("{},{},{},{:e}", i, j, k, f);
    }
}
