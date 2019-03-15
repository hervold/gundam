extern crate gundam;
#[macro_use]
extern crate log;
extern crate bio;
extern crate env_logger;
extern crate fishers_exact;
extern crate jobsteal;
extern crate ndarray;
extern crate suffix;

use std::collections::HashSet;
use bio::alphabets::dna::revcomp;

use gundam::*;
use bio::pattern_matching::pssm::{DNAMotif, Motif};
use env_logger::Builder as LogBuilder;
use fishers_exact::{fishers_exact, TestTails};
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
use suffix::SuffixTable;
use std::str;

const EPSILON: f32 = 1e-4;
const DELIMITER: u8 = b'$';

struct InputRec {
    idx1: usize,
    idx2: usize,
    gap: usize,
    motif_0: Vec<u8>,
    motif_1: Vec<u8>,
}

/// given motif, choose matching sequences and generate new motif representing mean
fn seqs_to_dyad<'a>(
    mut cpu_pool: &mut Pool,
    init: &DNAMotif,
    pos: &'a Vec<Vec<u8>>,
    neg: &'a Vec<Vec<u8>>,
) -> (f64, DyadMotif<'a, DNAMotif>) {
    let mut pos_v = init.eval_seqs(cpu_pool, pos);
    let mut neg_v = init.eval_seqs(cpu_pool, neg);
    let (pos_seqs, neg_ct, neg_seqs) = choose(&mut pos_v, &mut neg_v);

    info!("{} / {} seqs", pos_seqs.len(), pos.len());

    (scaled_fisher(pos_seqs.len(), pos_v.len(), neg_ct, neg_v.len()),
        DyadMotif {
        init: init.clone(),
        history: vec![MotifHistory::Init],
        motif: init.clone(),
        kmer_len: 0,
        gap_len: 0,
        pos_seqs: pos_seqs,
        neg_seqs: neg_seqs,
        score: f64::NAN,
    })
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
) -> (f32, f64, DNAMotif) {
    let (p_val, dyad) = seqs_to_dyad(&mut cpu_pool, init, pos, neg);
    let mean = dyad.refine_mean();

    (dist(&dyad.motif.scores, &mean.motif.scores), p_val, mean.motif)
}

fn mean_until_stable<'a>(
    mut cpu_pool: &mut Pool,
    init: &DNAMotif,
    pos: &'a Vec<Vec<u8>>,
    neg: &'a Vec<Vec<u8>>,
) -> (f64, DyadMotif<'a, DNAMotif>) {
    let (dist, p_val, motif) = em_cycle(&mut cpu_pool, init, pos, neg);
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
    let mut uniq: HashSet<Vec<u8>> = HashSet::new();

    let file = File::open(&args[1]).expect("can't open motifs file");
    for line in BufReader::new(&file).lines() {
        uniq.insert(line?.as_bytes().to_vec());
    }

    let mut motif_v: Vec<Vec<u8>> = uniq.into_iter().collect();
    motif_v.sort_by_key(|s| -1 * s.len() as isize);

    let mut not_subst_idx = vec![];
    let mut all_seqs: Vec<u8> = vec![];
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
        let motif = DNAMotif::from_degen(motif_v[idx].as_ref())?;
        let (p_val, dyad_final) = mean_until_stable(&mut pool, &motif, &pos, &neg);
        println!(
            "{},{},{:e}",
            String::from_utf8(motif.degenerate_consensus())?,
            String::from_utf8(dyad_final.motif.degenerate_consensus())?,
            p_val
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
