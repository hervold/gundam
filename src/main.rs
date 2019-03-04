extern crate gundam;
#[macro_use]
extern crate log;
extern crate bio;
extern crate chrono;
extern crate darwin_rs;
extern crate env_logger;
extern crate ndarray;

extern crate suffix;

use std::collections::HashSet;
use std::error::Error;

use bio::alignment::distance::hamming;
use bio::pattern_matching::pssm::Motif;
use bio::pattern_matching::pssm::{DNAMotif, ScoredPos};
use chrono::Local;
use darwin_rs::individual::Individual;
use env_logger::Builder as LogBuilder;
use gundam::*;
use ndarray::prelude::Array2;
use suffix::SuffixTable;

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::exit;
use std::str;

const DELIMITER: u8 = b'$';

fn extr_subseq<'a>(seq: &'a [u8], len: usize, sp: &ScoredPos) -> &'a [u8] {
    return &seq[sp.loc..sp.loc + len];
}

fn main() -> Result<(), Box<Error>> {
    //env_logger::init();
    let _ = LogBuilder::new()
        .parse(&env::var("RUST_LOG").unwrap_or_default())
        .init();

    let args = env::args().collect::<Vec<String>>();
    if args.len() != 4 {
        println!(
            "invalid # of args: {}.  usage: gundam <indices.txt> <pos.fa> <neg.fa>",
            args.len()
        );
        exit(1);
    }

    let file = File::open(&args[1]).expect("can't open index file");
    let idx_file = BufReader::new(&file);
    let indices: Vec<(usize, usize, usize, f64)> = idx_file
        .lines()
        .map(|line| {
            let a = line
                .as_ref()
                .expect("no line?")
                .split(",")
                .collect::<Vec<&str>>();
            (
                a[0].parse::<usize>().expect("first"),
                a[1].parse::<usize>().expect("second"),
                a[2].parse::<usize>().expect("third"),
                a[3].parse::<f64>().expect("fourth"),
            )
        })
        .collect();

    let mut uniq = HashSet::new();
    let all_dyads = DyadMotif::<DNAMotif>::motifs(indices, &args[2], &args[3], dyad::choose);
    for dyad in all_dyads.iter() {
        uniq.insert(dyad.refine_mean().motif.degenerate_consensus());
    }
    info!("{} unique motifs out of {}", uniq.len(), all_dyads.len());
    let mut motif_v: Vec<_> = uniq.iter().collect();
    motif_v.sort_by_key(|s| -1 * s.len() as isize);

    let mut all_seqs: Vec<u8> = vec![];
    let mut not_subst_idx = vec![];
    for (idx, motif_s) in motif_v.iter().enumerate() {
        let found = SuffixTable::new(str::from_utf8(all_seqs.as_ref())?)
            .contains(str::from_utf8(motif_s.as_ref())?);

        if !found {
            all_seqs.extend(motif_s.iter());
            all_seqs.push(DELIMITER);

            not_subst_idx.push(idx);
        }
    }

    for idx in not_subst_idx {
        println!("{}", str::from_utf8(motif_v[idx].as_ref())?);
    }

    Ok(())
}
