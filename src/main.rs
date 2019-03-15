extern crate gundam;
#[macro_use]
extern crate log;
extern crate bio;
extern crate chrono;
extern crate darwin_rs;
extern crate env_logger;
extern crate ndarray;

extern crate suffix;

use std::error::Error;

use bio::alignment::distance::hamming;
use bio::alphabets::dna::revcomp;
use bio::pattern_matching::pssm::{DNAMotif, Motif, ScoredPos};
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

    let pos = read_seqs(&args[2]);
    let neg = read_seqs(&args[3]);

    for idx in indices.into_iter() {
        // there should be only one ...
        for dyad in DyadMotif::<DNAMotif>::motifs(vec![idx], &pos, &neg, dyad::choose).iter() {
            let degen_before = dyad.motif.degenerate_consensus();
            let mean = dyad.refine_mean();
            let degen_after = mean.motif.degenerate_consensus();
            println!("{},{},{},{:e},{},{}",
                     idx.0,
                     idx.1,
                     idx.2,
                     idx.3,
                     str::from_utf8(degen_before.as_ref())?,
                     str::from_utf8(degen_after.as_ref())?,);
        }
    }

    /*
    let all_dyads = DyadMotif::<DNAMotif>::motifs(indices, &pos, &neg, dyad::choose);
    for dyad in all_dyads.iter() {
        uniq.insert(dyad.refine_mean().motif.degenerate_consensus());
    }
    info!("{} unique motifs out of {}", uniq.len(), all_dyads.len());
    let mut motif_v: Vec<Vec<u8>> = uniq.into_iter().collect();
    motif_v.sort_by_key(|s| -1 * s.len() as isize);

    let mut all_seqs: Vec<u8> = vec![];
    let mut not_subst_idx = vec![];
    for (idx, motif_s) in motif_v.iter().enumerate() {
        let found = {
            let suff_table = SuffixTable::new(str::from_utf8(all_seqs.as_ref())?);
            let s: &[u8] = motif_s.as_ref();
            suff_table.contains(str::from_utf8(s)?)
                || suff_table.contains(str::from_utf8(revcomp(s).as_slice())?)
        };
        if !found {
            all_seqs.extend(motif_s.iter());
            all_seqs.push(DELIMITER);

            not_subst_idx.push(idx);
        }
    }

    for idx in not_subst_idx {
        println!("{}", str::from_utf8(motif_v[idx].as_ref())?);
    }
     */
    Ok(())
}
