extern crate gundam;
#[macro_use]
extern crate log;
extern crate bio;
extern crate chrono;
extern crate env_logger;
extern crate ndarray;

extern crate suffix;

use std::error::Error;

use bio::alignment::distance::hamming;
use bio::alphabets::dna::revcomp;
use bio::pattern_matching::pssm::{DNAMotif, Motif, ScoredPos};
use chrono::Local;
//use darwin_rs::individual::Individual;
use env_logger::Builder as LogBuilder;
use gundam::*;
use ndarray::prelude::Array2;
use suffix::SuffixTable;

use std::cmp::{max, min};
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

    let b = b"AAAAGNNNNNNNNNNNNNNNNNNNNNAATC";
    let m = DNAMotif::from_degenerate(b).unwrap();
    info!(
        "{} - info: {}, thresh: {}",
        str::from_utf8(b).unwrap(),
        m.info_content(),
        passing_threshold(&m)
    );
    let b = b"TCTGAGAGGTGTTGAACGTTCTTGCTTCTG";
    let m = DNAMotif::from_degenerate(b).unwrap();
    info!(
        "{} - info: {}, thresh: {}",
        str::from_utf8(b).unwrap(),
        m.info_content(),
        passing_threshold(&m)
    );

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
    let indices: Vec<IndexRep> = idx_file
        .lines()
        .map(|line| {
            let a = line
                .as_ref()
                .expect("no line?")
                .split(",")
                .collect::<Vec<&str>>();

            IndexRep {
                kmer1: a[0].parse::<usize>().expect("1st"),
                gap: a[1].parse::<usize>().expect("2nd"),
                kmer2: a[2].parse::<usize>().expect("3rd"),
                kmer1_idx: a[3].parse::<usize>().expect("4th"),
                kmer2_idx: a[4].parse::<usize>().expect("5th"),
                p_val: a[5].parse::<f64>().expect("6th"),
            }
        })
        .collect();

    info!("finished reading file; found {} indices", indices.len());
    info!("~~~~ indices: {:#?}", &indices);

    let pos = read_seqs(&args[2]);
    let neg = read_seqs(&args[3]);

    // for each index, create a dyad from the index, and create a motif representing
    // the "mean" of matching sequences
    // then collect the degenerate representation of these motifs into a Vec
    let indices_len = indices.len();
    let log_every = max(1, min(10_000, indices_len / 50));
    info!("logging every {}", log_every);
    let mut degen_motifs = indices
        .into_iter()
        .enumerate()
        .filter_map(|(ct, idx)| {
            if ct % log_every == 0 {
                info!("-- processing index {} / {}", ct, indices_len);
            }

            let motif_v_from_idx =
                DyadMotif::<DNAMotif>::motifs(vec![idx], &pos, &neg, dyad::choose);
            match motif_v_from_idx.len() {
                0 => None,
                1 => {
                    if ct % log_every == 0 {
                        info!(
                            "-- motif - pos: {}/{}, neg: {}",
                            motif_v_from_idx[0].pos_seqs.len(),
                            pos.len(),
                            motif_v_from_idx[0].neg_seq_ct,
                        );
                    }
                    let mean = motif_v_from_idx[0].refine_mean();
                    if ct % log_every == 0 {
                        info!(
                            "-- index {}: {} -> {}",
                            ct,
                            str::from_utf8(
                                motif_v_from_idx[0].motif.degenerate_consensus().as_ref()
                            )
                            .unwrap(),
                            str::from_utf8(mean.motif.degenerate_consensus().as_ref()).unwrap()
                        );
                    }

                    Some(mean.motif.degenerate_consensus())
                }
                _ => unreachable!(),
            }
        })
        .collect::<Vec<_>>();

    // sort the Vec of degenerate motifs
    degen_motifs.sort();
    degen_motifs.sort_by_key(|s| -1 * s.len() as isize);

    info!("-- finished generating degenerate motifs");

    // use a suffix tree to find the longest version of each motif
    let mut not_subst_idx = vec![];
    let mut all_seqs: Vec<u8> = vec![];
    for (idx, motif_s) in degen_motifs.iter().enumerate() {
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

    info!("-- found uniqs");

    for i in not_subst_idx.into_iter() {
        println!("{}", str::from_utf8(degen_motifs[i].as_ref()).unwrap());
    }

    Ok(())
}
