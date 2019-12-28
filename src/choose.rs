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

    for motif_len in &[10, 15, 20, 25] {
        for rec in
            DyadMotif::<DNAMotif>::passing_kmers(*motif_len, args[1].as_str(), args[2].as_str())
        {
            println!(
                "{},{},{},{},{},{:e},{},{}",
                rec.kmer1,
                rec.gap,
                rec.kmer2,
                rec.kmer1_idx,
                rec.kmer2_idx,
                rec.p_val,
                Seq(GappedKmerCtr::<DNAMotif>::int_to_kmer(rec.kmer1, rec.kmer1_idx).as_ref()),
                Seq(GappedKmerCtr::<DNAMotif>::int_to_kmer(rec.kmer2, rec.kmer2_idx).as_ref()),
            );
        }
    }
}
