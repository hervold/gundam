extern crate gundam;
#[macro_use]
extern crate log;
extern crate bio;
extern crate env_logger;
extern crate jobsteal;
extern crate ndarray;

use bio::pattern_matching::pssm::DNAMotif;
use bio::pattern_matching::pssm::Motif;
use env_logger::Builder as LogBuilder;
use gundam::dyad::choose;
use gundam::*;
use jobsteal::{make_pool, BorrowSpliterator, IntoSpliterator, Pool, Spliterator};
use ndarray::prelude::{Array, Array2, AsArray};
use std::env;
use std::error::Error;
use std::f64;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::exit;

const EPSILON: f32 = 1e-4;

/// given motif, choose matching sequences and generate new motif representing mean
fn seqs_to_dyad<'a>(
    mut cpu_pool: &mut Pool,
    init: &DNAMotif,
    pos: &'a Vec<Vec<u8>>,
    neg: &'a Vec<Vec<u8>>,
) -> DyadMotif<'a, DNAMotif> {
    let mut pos_v = init.eval_seqs(cpu_pool, pos);
    let mut neg_v = init.eval_seqs(cpu_pool, neg);
    let (pos_seqs, neg_seqs) = choose(&mut pos_v, &mut neg_v).expect("motifs found bad one (1)");

    info!("{} / {} seqs", pos_seqs.len(), pos.len());

    DyadMotif {
        init: init.clone(),
        history: vec![MotifHistory::Init],
        motif: init.clone(),
        kmer_len: 0,
        gap_len: 0,
        pos_seqs: pos_seqs,
        neg_seqs: neg_seqs,
        score: f64::NAN,
    }
}

// FIXME: make this more generic ...
fn dist(a: &Array2<f32>, b: &Array2<f32>) -> f32 {
    a.iter()
        .zip(b.iter())
        .map(|(f, s)| {
            let v = f - s;
            v * v
        })
        .sum::<f32>()
        .sqrt()
}

// apply motif to find sequences; take mean; return new dyad + distance calculation
fn em_cycle<'a>(
    mut cpu_pool: &mut Pool,
    init: &DNAMotif,
    pos: &'a Vec<Vec<u8>>,
    neg: &'a Vec<Vec<u8>>,
) -> (f32, DNAMotif) {
    let dyad = seqs_to_dyad(&mut cpu_pool, init, pos, neg);
    let mean = dyad.refine_mean();

    (dist(&dyad.motif.scores, &mean.motif.scores), mean.motif)
}

fn mean_until_stable<'a>(
    mut cpu_pool: &mut Pool,
    init: &DNAMotif,
    pos: &'a Vec<Vec<u8>>,
    neg: &'a Vec<Vec<u8>>,
) -> DyadMotif<'a, DNAMotif> {
    let (dist, motif) = em_cycle(&mut cpu_pool, init, pos, neg);
    if dist < EPSILON {
        info!("found it!");
        seqs_to_dyad(&mut cpu_pool, &motif, pos, neg)
    } else {
        info!("recursing...");
        mean_until_stable(&mut cpu_pool, &motif, pos, neg)
    }
}

fn main() -> Result<(), Box<Error>> {
    let _ = LogBuilder::new()
        .parse(&env::var("RUST_LOG").unwrap_or_default())
        .try_init();

    warn!("starting..");

    let args = env::args().collect::<Vec<String>>();
    if args.len() != 4 {
        println!(
            "invalid # of args: {}.  usage: gundam_converge <motifs.txt> <pos.fa> <neg.fa>",
            args.len()
        );
        exit(1);
    }

    let pos = read_seqs(&args[2]);
    let neg = read_seqs(&args[3]);

    let mut pool = make_pool(*CPU_COUNT).unwrap();

    let file = File::open(&args[1]).expect("can't open motifs file");
    for line in BufReader::new(&file).lines() {
        let motif = DNAMotif::from_degen(line?.as_bytes())?;
        let dyad_final = mean_until_stable(&mut pool, &motif, &pos, &neg);
        println!(
            "{},{}",
            String::from_utf8(motif.degenerate_consensus())?,
            String::from_utf8(dyad_final.motif.degenerate_consensus())?
        );
        //let dyad_1 = dyad_0.refine_mean();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dist() {
        let a: Array2<f32> = Array2::zeros((3, 3));
        let b: Array2<f32> = Array2::from_elem((3, 3), 1.0);
        assert_eq!(dist(&a, &a), 0.0);
        assert!((3.0 - dist(&a, &b)).abs() < 0.0001);
    }
}
