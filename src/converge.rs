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
use std::collections::HashSet;

use bio::pattern_matching::pssm::{DNAMotif, Motif, PSSMError};
use env_logger::Builder as LogBuilder;
use fishers_exact::{fishers_exact, TestTails};
use jobsteal::{make_pool, BorrowSpliterator, IntoSpliterator, Pool, Spliterator};
use ndarray::prelude::{Array, Array2, AsArray};
use std::env;
use std::error::Error;
use std::f64;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::exit;
use std::str;
//use suffix::SuffixTable;

const DELIMITER: u8 = b'$';

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
    /*
    for line in BufReader::new(&file).lines() {
        uniq.insert(line?.as_bytes().to_vec());
    }*/

    let indices: Vec<(usize, usize, usize, f64)> = BufReader::new(file)
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

    let pos_idx = KmerIndex::new(&pos);
    let neg_idx = KmerIndex::new(&neg);

    /// stage 1: we have kmer pairs, which we'll
    /// (1) convert to DNAMotif's and use to pick sequnences, then
    /// (2) use those sequences to find a mean representation
    /// (3) uniq representations are checked for uniqueness
    /// ie, this is the first cycle of our EM
    let mut uniq: HashSet<Vec<u8>> = HashSet::new();
    for (i, j, k, _) in indices.into_iter() {
        let init: DNAMotif = (DyadMotif::<DNAMotif>::kmers_to_matrix(
            GappedKmerCtr::<DNAMotif>::int_to_kmer(KMER_LEN, i).as_slice(),
            k,
            GappedKmerCtr::<DNAMotif>::int_to_kmer(KMER_LEN, j).as_slice(),
        ))
        .into();

        let x: Result<(), Box<Error>> = (|| {
            let chosen_pos_idx = kmer_heur(&init, 0.99, &pos_idx)?;
            let chosen_pos = chosen_pos_idx
                .iter()
                .map(|i| pos[*i].clone())
                .collect::<Vec<_>>();
            let chosen_neg_idx = kmer_heur(&init, 0.99, &pos_idx)?;
            let chosen_neg = chosen_neg_idx
                .iter()
                .map(|i| neg[*i].clone())
                .collect::<Vec<_>>();

            let (p_val, dyad) = seqs_to_dyad(
                &mut pool,
                &init,
                &chosen_pos,
                pos.len(),
                &chosen_neg,
                neg.len(),
            );

            let degen_before = dyad.motif.degenerate_consensus();
            let mean = dyad.refine_mean();
            let degen_after = mean.motif.degenerate_consensus();

            uniq.insert(degen_after);
            Ok(())
        })();
    }
    info!("-- finished creating {} mean-motifs", uniq.len());

    let mut motif_v: Vec<Vec<u8>> = uniq.into_iter().collect();
    motif_v.sort();
    motif_v.sort_by_key(|s| -1 * s.len() as isize);

    let mut not_subst_idx = vec![];
    let mut all_seqs: Vec<u8> = vec![];
    for (idx, motif_s) in motif_v.iter().enumerate() {
        let found = {
            let suff_table = suffix::SuffixTable::new(str::from_utf8(all_seqs.as_ref())?);
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

    info!("-- after suffix tree, {} left", not_subst_idx.len());

    //let vals = vec![];
    let mut uniq: HashSet<Vec<u8>> = HashSet::new();
    for (num, idx) in not_subst_idx.iter().enumerate() {
        if num % 50 == 0 {
            info!("motif #{} / {}", num, not_subst_idx.len());
        }
        let motif_init = DNAMotif::from_degenerate(motif_v[*idx].as_ref())?;
        let mut final_pval: f64 = 0.0;
        let mut final_cts = (0, 0);
        let mut final_degen = vec![];

        match (|| -> Result<(), Box<Error>> {
            let mut chosen_pos = vec![];
            let mut chosen_neg = vec![];
            let mut motif = motif_init.clone();
            let mut p_val: f64 = 0.0;
            let mut hist = vec![];

            loop {
                let mut final_dyad = DyadMotif::<DNAMotif>::new();

                let chosen_pos_idx = kmer_heur(&motif, 0.99, &pos_idx)?;
                chosen_pos = chosen_pos_idx.iter().map(|i| pos[*i].clone()).collect();
                let chosen_neg_idx = kmer_heur(&motif, 0.99, &pos_idx)?;
                chosen_neg = chosen_neg_idx.iter().map(|i| neg[*i].clone()).collect();

                let t = mean_until_stable(
                    &mut pool,
                    &motif,
                    &chosen_pos,
                    pos.len(),
                    &chosen_neg,
                    neg.len(),
                );
                p_val = t.0;
                final_dyad = t.1;
                final_dyad.history = hist;

                if final_dyad.slide() {
                    motif = final_dyad.motif.clone();
                    hist = final_dyad.history.clone();
                } else {
                    final_pval = p_val;
                    final_degen = final_dyad.motif.degenerate_consensus();
                    final_cts = (final_dyad.pos_seqs.len(), final_dyad.neg_seq_ct);
                    break;
                }
            }

            Ok(())
        })() {
            Ok(()) => (),
            Err(_) => {
                final_pval = -1.0;
            }
        }
        /*
        vals.push((
            motif_init.degenerate_consensus(),
            final_degen,
            final_pval,
            final_cts.0,
            pos.len(),
            final_cts.1,
            neg.len(),
        ));*/

        if !uniq.contains(final_degen.as_slice()) {
            println!(
                "{},{},{:e},{},{},{},{}",
                String::from_utf8(motif_init.degenerate_consensus())?,
                str::from_utf8(final_degen.as_slice())?,
                final_pval,
                final_cts.0,
                pos.len(),
                final_cts.1,
                neg.len(),
            );
            uniq.insert(final_degen);
        }
    }

    Ok(())
}
