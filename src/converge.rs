extern crate bio;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate fishers_exact;
extern crate gundam;
extern crate itertools;
extern crate jobsteal;
extern crate ndarray;
extern crate suffix;
#[macro_use]
extern crate ergo_std;
use bio::alphabets::dna::revcomp;
use bio::pattern_matching::pssm::{DNAMotif, Motif, PSSMError};
use env_logger::Builder as LogBuilder;
use fishers_exact::{fishers_exact, TestTails};
use gundam::ctr::GappedKmerCtr;
use gundam::dyad::{choose, read_seqs, DyadMotif, MatrixPlus, MotifHistory};
use gundam::kmer_idx;
use gundam::kmer_idx::KmerIndex;
use gundam::*;
use gundam::{scaled_fisher, CPU_COUNT, KMER_LEN};
use itertools::join;
use itertools::Itertools;
use jobsteal::{make_pool, BorrowSpliterator, IntoSpliterator, Pool, Spliterator};
use ndarray::prelude::{Array, Array2, AsArray};
use std::cmp::min;
use std::collections::HashMap;
use std::collections::HashSet;
use std::env;
use std::error::Error;
use std::f64;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::process::exit;
use std::str;
use suffix::SuffixTable;

const KMER_SPLIT: usize = 5;
const DELIMITER: u8 = b'$';
const THRESH: f32 = 0.7;
const THRESH_CT: usize = 20;

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

    info!("pos.len(): {}, neg.len(): {}", pos.len(), neg.len());

    let pos_idx = KmerIndex::new(&pos);
    let neg_idx = KmerIndex::new(&neg);

    let mut pool = make_pool(*CPU_COUNT).unwrap();

    println!("idx,initial_motif,final_motif,init_info_content,threshold,values");
    for thresh_incr in 0..THRESH_CT {
        let threshold: f32 = THRESH + thresh_incr as f32 * (1.0 - THRESH) / THRESH_CT as f32;
        info!("~~ thresh: {}", threshold);

        let file = File::open(&args[1]).expect("can't open motifs file");
        let mut motif_init = DNAMotif::from_degenerate(b"A").unwrap();

        for (num, motif_s_) in BufReader::new(file).lines().enumerate() {
            let mut scores: HashMap<Vec<u8>, Vec<String>> = HashMap::new();

            let motif_s = motif_s_.unwrap();

            if num % 50 == 0 {
                info!("motif #{}: {}", num, motif_s.as_str());
            }
            motif_init = DNAMotif::from_degenerate(motif_s.as_bytes()).unwrap();
            let mut final_pval: f64 = 0.0;
            let mut final_cts = (0, 0);
            let mut final_degen = vec![];

            for window_num in 0..KMER_SPLIT {
                if num % 50 == 0 {
                    info!("~~ motif #{}, window #{}", num, window_num);
                }

                // we likely have different sequence counts, so the windows differ
                let pos_exclude = window_num * (pos.len() / KMER_SPLIT)
                    ..min((window_num + 1) * (pos.len() / KMER_SPLIT), pos.len());
                let neg_exclude = window_num * (neg.len() / KMER_SPLIT)
                    ..min((window_num + 1) * (neg.len() / KMER_SPLIT), neg.len());

                let kmer_pos = kmer_heur(&motif_init, &pos_idx, Some(&pos_exclude))?;
                let kmer_neg = kmer_heur(&motif_init, &neg_idx, Some(&neg_exclude))?;

                match (|| -> Result<(), Box<Error>> {
                    let mut motif = motif_init.clone();
                    let mut p_val: f64 = 0.0;
                    let mut hist = vec![];

                    // loop to handle "slide" operation
                    loop {
                        let mut final_dyad = DyadMotif::<DNAMotif>::new();
                        let mut pos_v = motif.eval_seqs(
                            &mut pool,
                            pos.iter().enumerate().filter_map(|(i, s)| {
                                if kmer_pos.contains(&i) {
                                    Some(s.as_ref())
                                } else {
                                    None
                                }
                            }),
                        );
                        let mut neg_v = motif.eval_seqs(
                            &mut pool,
                            neg.iter().enumerate().filter_map(|(i, s)| {
                                if kmer_neg.contains(&i) {
                                    Some(s.as_ref())
                                } else {
                                    None
                                }
                            }),
                        );

                        let (chosen_pos, _, chosen_neg) = choose(threshold, &mut pos_v, &mut neg_v);
                        let t = mean_until_stable(
                            &mut pool,
                            &motif,
                            &chosen_pos,
                            pos.len(),
                            &chosen_neg,
                            neg.len(),
                            threshold,
                        );
                        p_val = t.0;
                        final_dyad = t.1;
                        final_dyad.history = hist;

                        info!(
                            "~~ after stable: {}",
                            Seq(final_dyad.motif.degenerate_consensus().as_ref())
                        );

                        if final_dyad.slide() {
                            motif = final_dyad.motif.clone();
                            hist = final_dyad.history.clone();
                        } else {
                            final_pval = p_val;
                            final_degen = final_dyad.motif.degenerate_consensus();
                            final_cts = (final_dyad.pos_seqs.len(), final_dyad.neg_seq_ct);
                            info!(
                                "~~ hist for {}: {:?}",
                                Seq(final_degen.as_ref()),
                                &final_dyad.history
                            );
                            break;
                        }
                    }

                    Ok(())
                })() {
                    Ok(()) => (),
                    Err(e) => {
                        info!("~~~~~ ERROR: {:?}", e);
                        final_pval = -1.0;
                    }
                }
                info!(
                    "~~ #{}, window {} - final degen: {}",
                    num,
                    window_num,
                    Seq(final_degen.as_ref())
                );

                // we should have a motif for this batch of sequences, excluding the window - ie, most
                // sequences were used to train
                // next calculate its score using the excluded sequences, ie, the validation set
                if final_pval != -1.0 {
                    let motif = DNAMotif::from_degenerate(final_degen.as_ref()).unwrap();
                    let (pval, _) = seqs_to_dyad(
                        &mut pool,
                        &motif,
                        SeqRefList::Without(&mut pos_exclude.map(|i| pos[i].as_ref())),
                        pos.len(),
                        SeqRefList::Without(&mut neg_exclude.map(|i| neg[i].as_ref())),
                        neg.len(),
                        Some(THRESH), //passing_threshold(&motif)),
                    );
                    scores
                        .entry(motif.degenerate_consensus())
                        .or_insert_with(|| vec![])
                        .push(format!("motif-{}|window-{}|{:.5e}", num, window_num, pval));
                }
            }
            for (degen, pvals) in scores {
                println!(
                    "{},{},{},{},{},{}",
                    num,
                    motif_s,
                    Seq(degen.as_ref()),
                    motif_init.info_content(),
                    threshold,
                    pvals.iter().join(",")
                );
            }
        }
    }
    Ok(())
}
