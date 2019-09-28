extern crate bio;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate fishers_exact;
extern crate gundam;
extern crate jobsteal;
extern crate ndarray;
extern crate suffix;
#[macro_use]
extern crate ergo_std;

use gundam::ctr::GappedKmerCtr;
use gundam::dyad::{choose, read_seqs, DyadMotif, MatrixPlus, MotifHistory};
use gundam::kmer_idx;
use gundam::kmer_idx::KmerIndex;
use gundam::*;
use gundam::{scaled_fisher, CPU_COUNT, KMER_LEN};

use bio::alphabets::dna::revcomp;
use bio::io::fasta;
use bio::pattern_matching::pssm::{DNAMotif, Motif, PSSMError, ScoredPos};
use env_logger::Builder as LogBuilder;
use fishers_exact::{fishers_exact, TestTails};
use jobsteal::{
    make_pool, BorrowSpliterator, BorrowSpliteratorMut, IntoSpliterator, Pool, Spliterator,
};
use ndarray::prelude::{Array, Array2, AsArray};
use std::collections::HashSet;
use std::env;
use std::error::Error;
use std::f64;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::mem;
use std::process::exit;
use std::str;

fn tally_hits(
    pool: &mut Pool,
    motifs: &Vec<(Vec<u8>, DNAMotif)>,
    fname: &str,
) -> (Vec<usize>, usize) {
    let mut seqs = fasta::Reader::from_file(fname)
        .expect(format!("couldn't open {}", fname).as_str())
        .records()
        .enumerate()
        .map(|(idx, rec_)| (idx, ScoredPos::default(), rec_, 0))
        .collect::<Vec<_>>();

    let seq_ct = seqs.len();
    seqs.split_iter_mut().for_each(
        &pool.spawner(),
        |&mut (seq_idx, ref mut sp, ref rec_, ref mut best_idx)| {
            if seq_idx % 200 == 0 {
                info!("#{} / {}", seq_idx, seq_ct);
            }
            for (idx, (_, m)) in motifs.iter().enumerate() {
                let rec: &fasta::Record = rec_.as_ref().expect("rec?");
                let mut next_score: ScoredPos = m.score(rec.seq()).expect("score?");
                if next_score.sum > sp.sum {
                    let _ = mem::replace(sp, next_score);
                    let _ = mem::replace(best_idx, idx);
                }
            }
        },
    );

    let mut cts = vec![0_usize; motifs.len()];
    for (_, _, _, best_idx) in seqs {
        cts[best_idx] += 1;
    }

    (cts, seq_ct)
}

fn main() -> Result<(), Box<Error>> {
    let _ = LogBuilder::new()
        .parse(&env::var("RUST_LOG").unwrap_or_default())
        .try_init();

    warn!("starting..");

    let args = env::args().collect::<Vec<String>>();
    if args.len() < 3 || args.len() > 4 {
        println!(
            "invalid # of args: {}.  usage: gundam_match <motifs.txt> <queries.fa> [<neg_queries.fa>]",
            args.len()
        );
        exit(1);
    }

    let mut motifs = vec![];
    for line in BufReader::new(&File::open(&args[1])?).lines() {
        let degen = line?.into_bytes();
        let m = DNAMotif::from_degenerate(&degen[..])?;
        motifs.push((degen, m));
    }

    let mut pool = make_pool(*CPU_COUNT).unwrap();

    let pos = tally_hits(&mut pool, &motifs, args[2].as_str());
    let neg = if args.len() > 3 {
        Some(tally_hits(&mut pool, &motifs, args[3].as_str()))
    } else {
        None
    };

    for (idx, (degen, _)) in motifs.iter().enumerate() {
        match neg {
            Some(ref n) => {
                let pval = scaled_fisher(pos.0[idx], n.0[idx], pos.1, n.1);
                println!(
                    "{},{},{},{},{},{:e}",
                    str::from_utf8(degen)?,
                    pos.0[idx],
                    pos.1,
                    n.0[idx],
                    n.1,
                    pval
                );
            }
            None => {
                println!("{},{},{}", str::from_utf8(degen)?, pos.0[idx], pos.1);
            }
        }
    }

    Ok(())
}
