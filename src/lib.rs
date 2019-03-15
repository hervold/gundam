#[macro_use(s)]
extern crate ndarray;
extern crate bio;
extern crate darwin_rs;
extern crate jobsteal;
extern crate rand;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate fishers_exact;
extern crate num_cpus;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate serde_derive;

extern crate serde;
extern crate serde_json;
use fishers_exact::{fishers_exact, TestTails};
use std::cmp::{max, min};

use std::ffi::{CStr, CString};
use std::fmt;
use std::mem;
use std::os::raw::{c_char, c_void};

use darwin_rs::{Individual, Population, PopulationBuilder, SimulationBuilder};
use jobsteal::{make_pool, BorrowSpliteratorMut, Pool, Spliterator};
use std::f64;
use std::str;
//use darwin_rs::select::MaximizeSelector;

use bio::io::fasta;
use bio::pattern_matching::pssm::{DNAMotif, Motif, ScoredPos};
use ndarray::prelude::{Array, Array2, AsArray};
use rand::Rng;

pub mod ctr;
use ctr::*;

pub mod dyad;
pub use dyad::*;

//pub use dyad::find_motifs;

pub const KMER_LEN: usize = 5;
const MIN_GAP: usize = 0;
const MAX_GAP: usize = 20;
const MUT_INCR: f32 = 0.2;
const MIN_SCORE: f32 = 0.9;

lazy_static! {
    pub static ref CPU_COUNT: usize = num_cpus::get();
}

type ScoredSeqs = Vec<(String, f32)>;

/// P-values returned by the Fisher exact test don't change meaningfully as the values get larger.
/// eg, both [100, 200, 10_000, 10_000] and [1_000, 2_000, 100_000, 100_000] yield a P-value well
/// below our cutoff.  therefore, we can safely scale the values down if they're above some arbitrary
/// threshold.
pub fn scaled_fisher(_ct1: usize, _tot1: usize, _ct2: usize, _tot2: usize) -> f64 {
    let (ct1, tot1) = if _tot1 as f64 <= 1e4 {
        (_ct1 as i32, _tot1 as i32)
    } else {
        let scale = 1e4 / _tot1 as f64;
        (max(1, (scale * _ct1 as f64) as i32), 10_000)
    };

    let (ct2, tot2) = if _tot2 as f64 <= 1e4 {
        (_ct2 as i32, _tot2 as i32)
    } else {
        let scale = 1e4 / _tot2 as f64;
        (max(1, (scale * _ct2 as f64) as i32), 10_000)
    };

    fishers_exact(&[ct1, ct2, tot1, tot2], TestTails::One)
}

/*

interface ideas:

pub fn export_scores(dyad: &DyadMotif) -> (ScoredSeqs, ScoredSeqs)

pub fn get_dyad(*Vec<DyadMotif>, idx) -> DyadMotif (serde)

pub fn refine_GA(*Vec..., idx) -> idx of new

pub fn refine_mean(*Vec..., idx) -> idx of new

pub fn winnow_seqs(*Vec, idx, Vec<idx>)

pub fn shannon_entropy( Vec<seq> ) -> f32
*/

/*
#[no_mangle]
pub extern "C" fn release_str(somestr: *mut c_char) {
    unsafe {
        CString::from_raw(somestr);
    }
}

/// wrapper for DyadMotif::motifs
#[no_mangle]
pub extern "C" fn read_kmers(
    _kmers_s: *const c_char,
    _pos_fname: *const c_char,
    _neg_fname: *const c_char,
) -> *mut c_void {
    let kmers_s = unsafe { CStr::from_ptr(_kmers_s).to_string_lossy().into_owned() };
    let pos_fname = unsafe { CStr::from_ptr(_pos_fname).to_string_lossy().into_owned() };
    let neg_fname = unsafe { CStr::from_ptr(_neg_fname).to_string_lossy().into_owned() };
    let kmers: Vec<(usize, usize, usize, f64)> =
        serde_json::from_str(&kmers_s).expect("deser kmers");
    let motifs: Box<Vec<DyadMotif<DNAMotif>>> =
        Box::new(DyadMotif::<DNAMotif>::motifs(kmers, &pos_fname, &neg_fname, choose));
    Box::into_raw(motifs) as *mut c_void
}


/// input:
///    _dyads - Vec<DyadMotif> as returned by read_kmers
///    idx - u32
/// output:
///    DyadMotif at position idx, encoded as JSON
#[no_mangle]
pub extern "C" fn get_dyad(_dyads: *const c_void, idx: u32) -> *const c_char {
    let dyads: &Vec<DyadMotif<DNAMotif>> = unsafe { mem::transmute(_dyads) };

    CString::new(serde_json::to_string(&dyads[idx as usize]).expect(
        "get_dyad - ser",
    )).expect("get_dyad")
        .into_raw()
}

/// input:
///    _dyads - Vec<DyadMotif> as returned by read_kmers
/// output:
///    length of vec
#[no_mangle]
pub extern "C" fn get_len(_dyads: *const c_void) -> u32 {
    let dyads: &Vec<DyadMotif<DNAMotif>> = unsafe { mem::transmute(_dyads) };
    dyads.len() as u32
}


/// input:
///    _dyads - Vec<DyadMotif> as returned by read_kmers
///    idx - u32
/// output:
///    information content of motif
#[no_mangle]
pub extern "C" fn info_content(_dyads: *const c_void, idx: u32) -> f32 {
    let dyads: &Vec<DyadMotif<DNAMotif>> = unsafe { mem::transmute(_dyads) };

    dyads[idx as usize].motif.info_content()
}


/// input:
///    _dyads - Vec<DyadMotif> as returned by read_kmers
///    idx - u32
/// output:
///    information content of motif
#[no_mangle]
pub extern "C" fn show_motif(_dyads: *const c_void, idx: u32) -> *const c_char {
    let dyads: &Vec<DyadMotif<DNAMotif>> = unsafe { mem::transmute(_dyads) };

    CString::new(dyads[idx as usize].show_motif())
        .expect("get_dyad")
        .into_raw()
}


/// creates a new DyadMotif and appends it to the vec, returning
/// indes of new motif
/// input:
///    _dyads - Vec<DyadMotif> as returned by read_kmers
///    idx - u32
/// output:
///    index of new motif
#[no_mangle]
pub extern "C" fn simple_mean(_dyads: *const c_void, idx: u32) -> u32 {
    let dyads: &mut Vec<DyadMotif<DNAMotif>> = unsafe { mem::transmute(_dyads) };

    let new = dyads[idx as usize].refine_mean();
    let new_idx = dyads.len();
    dyads.push(new);
    new_idx as u32
}

*/

pub struct Seq<'a>(pub &'a [u8]);
impl<'a> fmt::Display for Seq<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", str::from_utf8(self.0).unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fishers_exact::{fishers_exact, TestTails};
    const MOTIF: &'static [u8] = b"GGCCTAGCCATG";
    const POS_FNAME: &'static str = "pos.fa";
    const NEG_FNAME: &'static str = "neg.fa";

    #[test]
    #[ignore]
    fn kmers_to_m() {
        let m = DyadMotif::<DNAMotif>::kmers_to_matrix(b"ATGC", 1, b"ATGC");
        let expected = Array::from_vec(vec![
            0.8, 0.05, 0.05, 0.05, 0.05, 0.8, 0.05, 0.05, 0.05, 0.05, 0.8, 0.05, 0.05, 0.05, 0.05,
            0.8, 0.25, 0.25, 0.25, 0.25, 0.8, 0.05, 0.05, 0.05, 0.05, 0.8, 0.05, 0.05, 0.05, 0.05,
            0.8, 0.05, 0.05, 0.05, 0.05, 0.8,
        ])
        .into_shape((9, 4))
        .unwrap();
        println!("diff: {:?}", m.clone() - expected.clone());
        assert_eq!(m, expected);
    }

    #[test]
    #[ignore]
    fn test_one() {
        let motif = DNAMotif::from(DyadMotif::<DNAMotif>::kmers_to_matrix(
            b"ATAGG", MAX_GAP, b"CCATG",
        ));
        println!("score for present: {:?}", motif.score(b"GGAACGAAGTCCGTAGGGTCCATAGGAAAACCACTATGGGGCAGGATAATCATTAAAGGTCACTCGGTCGAGGCACAGATTGTGAGGAAGATGTAGGGGACCGTCGTTAAACCTAACGGACGGCTACACGGTTGTTGAAATGTCCCCCCCTTTTGCATTTTTCCTATGGGCGGCGACATAAAACTCGCAGACGAAGTTGGATATCTCCCGAATACGTGGACCGGCAGCATAACCAGACAAACGGGTAACTAACGTATGAGTGTGTCCAGCCACCATCCATAGGAAGTCCCATGAGTGAGCTTGATGATGTGAGGGCATGACATGTGCGGAAAACGAAGAACTAGGACCATAATGCAGGGCGACCTGCGCTCGAAACTCTGGATTACCATTTCCGCGGCCTAATATGGATCTCCTGTGTCTCGGATCCTTCAGGTCGACGTTCGGATCATACATGGGACTACAACGTGTCGATAGACCGCCAGACCTACACAAAGCATGCA".iter()));
    }

    #[test]
    #[ignore]
    fn test_find() {
        let v = DyadMotif::<DNAMotif>::passing_kmers(POS_FNAME, NEG_FNAME);
        let pos = read_seqs(POS_FNAME);
        let neg = read_seqs(NEG_FNAME);
        let dyads: Vec<DyadMotif<DNAMotif>> = DyadMotif::<DNAMotif>::motifs(v, &pos, &neg, choose);
        //let new_dyad = dyads[0].refine(100);
    }

    #[test]
    #[ignore]
    fn print_kmers() {
        for i in 0..MOTIF.len() - KMER_LEN {
            println!(
                "@@ from motif, kmer {} -> {}",
                str::from_utf8(&MOTIF[i..i + KMER_LEN]).unwrap(),
                GappedKmerCtr::<DNAMotif>::kmer_to_int(&MOTIF[i..i + KMER_LEN])
            );
        }
    }
    #[test]
    fn fisher() {
        println!(
            "fisher: {:?}",
            fishers_exact(&[100, 200, 5000, 10000], TestTails::One)
        );
    }
}
