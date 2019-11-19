//use darwin_rs::{Individual, Population, PopulationBuilder, SimulationBuilder};
use bio::alignment::distance::hamming;
use bio::io::fasta;
use bio::pattern_matching::pssm::{Motif, PSSMError, ScoredPos};
use jobsteal::{make_pool, BorrowSpliterator, IntoSpliterator, Pool, Spliterator};
use ndarray::prelude::{Array, Array2};
use rand;
use rand::Rng;
use std::cmp::{max, min};
use std::ops::Range;
use std::str;
use std::{f64, usize};

use super::*;
use ctr::*;

const P_CUTOFF: f64 = 0.001;
const MODAL_MAX_SEQS: usize = 400;

#[derive(Debug, Clone, PartialEq)]
pub enum MotifHistory {
    Init,
    Mutate,
    Mean,
    MeanDiff,
    Mode,
    Reset,
    ExpandLeft,
    ExpandRight,
    ExpandBoth,
    ContractLeft,
    ContractRight,
    ContractBoth,
    Final,
}

#[derive(Debug, Clone)]
pub struct DyadMotif<'a, M>
where
    M: Motif + Clone,
{
    /// initial state based on kmers
    pub init: M,
    /// history of actions that produced this motif
    pub history: Vec<MotifHistory>,
    /// weights updated by GA
    pub motif: M,
    /// kmer lens
    pub kmer_len: (usize, usize),
    /// gap len
    pub gap_len: usize,
    /// sequences matching our motif
    pub pos_seqs: Vec<(&'a [u8], ScoredPos)>,
    /// sequences representing background
    pub neg_seqs: Vec<(&'a [u8], ScoredPos)>,
    /// number of passing negative sequences
    pub neg_seq_ct: usize,
    /// score - sum(neg) / sum(pos)
    pub score: f64,
}

impl<'a, M> DyadMotif<'a, M>
where
    M: 'a + Motif + Clone + Sync + Send + From<Array2<f32>>,
{
    pub fn new() -> DyadMotif<'static, DNAMotif> {
        let m = DNAMotif::from_degenerate(b"A").unwrap();
        DyadMotif {
            init: m.clone(),
            history: vec![],
            motif: m.clone(),
            kmer_len: (0, 0),
            gap_len: 0,
            pos_seqs: vec![],
            neg_seqs: vec![],
            neg_seq_ct: 0,
            score: 0.0,
        }
    }

    fn fasta_to_ctr_p(
        fname: &str,
        kmer_len: usize,
        min_gap: usize,
        max_gap: usize,
    ) -> (GappedKmerCtr<M>, usize) {
        let mut ctr = GappedKmerCtr::new(kmer_len, min_gap, max_gap);
        let mut tot = 0;

        for _rec in fasta::Reader::from_file(fname)
            .expect(format!("trouble opening {}", fname).as_str())
            .records()
        {
            let rec = _rec.expect("couldn't unwrap record");
            ctr.update_with_seq(rec.seq());
            tot += 1;
        }

        (ctr, tot)
    }

    fn fasta_to_ctr(fname: &str) -> (GappedKmerCtr<M>, usize) {
        DyadMotif::fasta_to_ctr_p(fname, KMER_LEN, MIN_GAP, MAX_GAP)
    }

    pub fn kmers_to_matrix(kmer1: &[u8], gap_len: usize, kmer2: &[u8]) -> Array2<f32> {
        let mut m = Array2::from_elem((kmer1.len() + gap_len + kmer2.len(), 4), 0.05);
        for i in 0..kmer1.len() {
            m[[i, M::lookup(kmer1[i]).expect("kmers_to_matrix")]] = 0.85;
        }
        // set gap to N, ie, equal weights
        for i in 0..gap_len + 1 {
            for j in 0..4 {
                m[[kmer1.len() + i, j]] = 0.25;
            }
        }
        for i in 0..kmer2.len() {
            m[[
                kmer1.len() + gap_len + i,
                M::lookup(kmer2[i]).expect("kmers_to_matrix"),
            ]] = 0.85;
        }
        m
    }

    /// generate kmers, tablulate, and apply Fisher exact test
    pub fn passing_kmers(
        pos_fname: &str,
        neg_fname: &str,
        specs: &Vec<(usize, usize, usize)>,
    ) -> Vec<(usize, usize, usize, f64)> {
        let cutoff = P_CUTOFF / 2.0_f64.powi(2 * KMER_LEN as i32);
        info!("-- using cutoff: {:e}", cutoff);

        let mut pool = make_pool(*CPU_COUNT).unwrap();
        debug!("-- passing_kmers - created pool of {} threads", *CPU_COUNT);
        let mut passing = Vec::new();

        let default = vec![(KMER_LEN, MIN_GAP, MAX_GAP)];
        for &(kmer_len, min_gap, max_gap) in (if specs.len() == 0 { &default } else { specs }) {
            let (pos, pos_ct) =
                DyadMotif::<DNAMotif>::fasta_to_ctr_p(pos_fname, kmer_len, min_gap, max_gap);
            let (neg, neg_ct) =
                DyadMotif::<DNAMotif>::fasta_to_ctr_p(neg_fname, kmer_len, min_gap, max_gap);
            let (width, height, gap) = pos.ctr.dim();

            let mut indices: Vec<(usize, Vec<(usize, usize, f64)>)> =
                (0..width).map(|i| (i, vec![])).collect();
            indices
                .split_iter_mut()
                .for_each(&pool.spawner(), |&mut (ref i, ref mut v)| {
                    //for i in 0..width {
                    for j in 0..height {
                        for k in 0..gap {
                            if pos.ctr[[*i, j, k]] > neg.ctr[[*i, j, k]] {
                                let p = scaled_fisher(
                                    pos.ctr[[*i, j, k]],
                                    pos_ct,
                                    neg.ctr[[*i, j, k]],
                                    neg_ct,
                                );
                                if p < P_CUTOFF {
                                    v.push((j, k, p));
                                };
                            }
                        }
                    }
                });
            for (i, v) in indices {
                passing.extend(v.iter().map(|&(j, k, p)| (i, j, k, p)));
            }
        }
        passing
    }

    pub fn motifs<F>(
        chosen: Vec<(usize, usize, usize, f64)>,
        pos: &'a Vec<Vec<u8>>,
        neg: &'a Vec<Vec<u8>>,
        chooser: F,
    ) -> Vec<DyadMotif<'a, DNAMotif>>
    where
        F: Fn(
            f32,
            &mut Vec<(&'a [u8], ScoredPos)>,
            &mut Vec<(&'a [u8], ScoredPos)>,
        ) -> (
            Vec<(&'a [u8], ScoredPos)>,
            usize,
            Vec<(&'a [u8], ScoredPos)>,
        ),
    {
        debug!("using {} cpus", *CPU_COUNT);
        let mut pool = make_pool(*CPU_COUNT).unwrap();
        let mut dyads = Vec::new();
        for (idx, &(i, j, k, _)) in chosen.iter().enumerate() {
            /*
            if idx % 500 == 0 {
                info!("creating dyad #{} / {}", idx, chosen.len());
            }*/

            let init: DNAMotif = (DyadMotif::<DNAMotif>::kmers_to_matrix(
                GappedKmerCtr::<DNAMotif>::int_to_kmer(KMER_LEN, i).as_slice(),
                k,
                GappedKmerCtr::<DNAMotif>::int_to_kmer(KMER_LEN, j).as_slice(),
            ))
            .into();

            // drop invalid motifs
            if init.min_score == init.max_score {
                debug!(
                    "skipping motif: {}",
                    String::from_utf8(init.degenerate_consensus()).expect("show_motif Q")
                );
                continue;
            } else {
                debug!(
                    "good motif {}<{}>{}: {}",
                    i,
                    k,
                    j,
                    String::from_utf8(init.degenerate_consensus()).expect("show_motif R")
                );
            }

            let mut pos_v = init.eval_seqs(&mut pool, pos.iter().map(|s| s.as_ref()));
            let mut neg_v = init.eval_seqs(&mut pool, neg.iter().map(|s| s.as_ref()));
            let (pos_seqs, neg_ct, neg_seqs) =
                chooser(passing_threshold(&init), &mut pos_v, &mut neg_v);

            debug!(
                "DyadMotif::motifs - {} / {} seqs used",
                pos_seqs.len(),
                pos_v.len()
            );

            let copy = init.clone();

            let mut d = DyadMotif {
                init: init,
                history: vec![MotifHistory::Init],
                motif: copy,
                kmer_len: (KMER_LEN, KMER_LEN),
                gap_len: MIN_GAP + k,
                pos_seqs: pos_seqs,
                neg_seqs: neg_seqs,
                neg_seq_ct: neg_ct,
                score: f64::NAN,
            };
            //d.calculate_fitness();
            dyads.push(d);
        }
        dyads
    }

    /// simple means-based refinement
    pub fn refine_mean(&self) -> DyadMotif<M> {
        let len = self.motif.len();
        let new_scores = tally(&self.motif, self.pos_seqs.iter().map(|(ref seq, _)| *seq));

        let m: M = new_scores.into();

        let mut hist = self.history.to_owned();
        hist.push(MotifHistory::Mean);
        let mut d: DyadMotif<M> = DyadMotif {
            init: self.init.clone(),
            history: hist,
            motif: m,
            kmer_len: self.kmer_len,
            gap_len: self.gap_len,
            pos_seqs: self.pos_seqs.to_vec(),
            neg_seqs: self.neg_seqs.to_vec(),
            neg_seq_ct: self.neg_seq_ct,
            score: f64::NAN,
        };

        //d.calculate_fitness();
        d
    }

    /// tally base-matches for positive seqs as in refine_mean, but also tally
    /// negetive seqs and divide positive tallies by negative before normalizing
    pub fn refine_meandiv(&self) -> DyadMotif<M> {
        let pos_tally = tally(&self.motif, self.pos_seqs.iter().map(|(ref seq, _)| *seq));
        let neg_tally = tally(&self.motif, self.neg_seqs.iter().map(|(ref seq, _)| *seq));

        let m: M = (pos_tally / neg_tally).into();
        //m.normalize_scores();
        //m.calc_minmax();

        info!(
            "meandiv: {:?}",
            String::from_utf8(m.degenerate_consensus()).unwrap()
        );

        let mut hist = self.history.to_owned();
        hist.push(MotifHistory::MeanDiff);
        let mut d = DyadMotif {
            init: self.init.clone(),
            history: hist,
            motif: m,
            kmer_len: self.kmer_len,
            gap_len: self.gap_len,
            pos_seqs: self.pos_seqs.to_vec(),
            neg_seqs: self.neg_seqs.to_vec(),
            neg_seq_ct: self.neg_seq_ct,
            score: f64::NAN,
        };

        //d.calculate_fitness();
        d
    }

    /// do a simple all-against-all Hamming distance comparison of pos_seqs,
    /// choose a sequence to represent a larger cluster (the "mode"), and cross
    /// the PSSM with that "modal" sequence
    /// only uses the first MODAL_MAX_SEQS sequences from self.pos_seqs, so any sorting of
    /// self.pos_seqs will bias outcome
    pub fn refine_via_mode(&self) -> DyadMotif<M> {
        let pwm_len = self.motif.len();
        let seq_ct = min(self.pos_seqs.len(), MODAL_MAX_SEQS);
        let mut diffs: Array2<usize> = Array2::zeros((seq_ct, seq_ct));
        for i in 0..seq_ct {
            for j in i..seq_ct {
                let (ref a_seq, ScoredPos { ref loc, .. }) = self.pos_seqs[i];
                let a = &a_seq[*loc..*loc + pwm_len];
                let (ref b_seq, ScoredPos { ref loc, .. }) = self.pos_seqs[j];
                let b = &b_seq[*loc..*loc + pwm_len];

                diffs[[i, j]] = hamming(a, b) as usize;
            }
        }
        let mut best_i = 0;
        let mut best_score = usize::MAX;
        for i in 0..seq_ct {
            let mut dists: Vec<usize> = (0..seq_ct)
                .map(|j| diffs[if i <= j { [i, j] } else { [j, i] }])
                .collect();
            dists.sort();
            let score: usize = dists[0..seq_ct / 4].iter().sum();
            if score < best_score {
                best_i = i;
                best_score = score;
            }
        }
        let (ref seq, ScoredPos { ref loc, .. }) = self.pos_seqs[best_i];

        let mut m = self.motif.get_scores().clone();
        for (i, b) in seq[*loc..*loc + pwm_len].iter().enumerate() {
            m[[i, M::lookup(*b).expect("refine_via_mode")]] += 0.5;
        }
        let motif: M = m.into();
        info!(
            "mode: {:?}",
            String::from_utf8(motif.degenerate_consensus()).unwrap()
        );

        let mut pos: Vec<(&[u8], ScoredPos)> = self
            .pos_seqs
            .iter()
            .map(|&(ref seq, _)| (seq.clone(), motif.score(*seq).expect("mode pos")))
            .collect();
        pos.sort_by(|&(_, ref score_a), &(_, ref score_b)| {
            score_b.sum.partial_cmp(&score_a.sum).expect("float sort")
        });
        let mut cutoff = 0;
        for (i, &(_, ref score)) in pos.iter().enumerate() {
            if score.sum <= 0.9 {
                break;
            }
            cutoff = i;
        }

        info!("mode cutoff of {} on {}", cutoff, pos.len());

        let mut neg: Vec<(&[u8], ScoredPos)> = self
            .neg_seqs
            .iter()
            .map(|&(ref seq, _)| (seq.clone(), motif.score(*seq).expect("mode neg")))
            .collect();
        neg.sort_by(|&(_, ref score_a), &(_, ref score_b)| {
            score_b.sum.partial_cmp(&score_a.sum).expect("float sort")
        });
        let mut hist = self.history.clone();
        hist.push(MotifHistory::Mode);
        let mut d = DyadMotif {
            init: self.motif.clone(),
            history: hist,
            motif: motif,
            kmer_len: self.kmer_len,
            gap_len: self.gap_len,
            pos_seqs: pos.iter().take(cutoff).cloned().collect(),
            neg_seqs: neg.iter().take(cutoff).cloned().collect(),
            neg_seq_ct: self.neg_seq_ct,
            score: f64::NAN,
        };
        //d.calculate_fitness();
        d
    }

    /*
    pub fn refine(&self, mut_ct: usize) -> DyadMotif<M> {
        self.refine_GA(mut_ct)
    }*/

    /// stringified self.motif.degenerate_consensus()
    pub fn show_motif(&self) -> String {
        String::from_utf8(self.motif.degenerate_consensus()).expect("show_motif")
    }

    /// expand or contract a dyad based upon degenerate representation:
    /// 1. if either 3' or 5' end is N, "cap" it: truncate and add a ContractX entry to the history vec
    /// 2. if not, lengthen by 1 N, and return for refinement
    pub fn slide(&mut self) -> bool {
        if let Some(entry) = self.history.last() {
            if *entry == MotifHistory::Final {
                return false;
            }
        }

        let left_term = match self
            .history
            .iter()
            .find(|e| **e == MotifHistory::ContractLeft)
        {
            Some(_) => true,
            None => false,
        };
        let right_term = match self
            .history
            .iter()
            .find(|e| **e == MotifHistory::ContractRight)
        {
            Some(_) => true,
            None => false,
        };

        if left_term && right_term {
            self.history.push(MotifHistory::Final);
            return false;
        }

        let mut new_len = self.motif.len();
        let repr = self.motif.degenerate_consensus();
        let left_incr = if left_term {
            0
        } else if repr[0] == b'N' {
            self.history.push(MotifHistory::ContractLeft);
            -1
        } else {
            self.history.push(MotifHistory::ExpandLeft);
            1
        };

        let right_incr = if right_term {
            0
        } else if let Some(b) = repr.last() {
            if *b == b'N' {
                self.history.push(MotifHistory::ContractRight);
                -1
            } else {
                self.history.push(MotifHistory::ExpandRight);
                1
            }
        } else {
            unreachable!();
        };

        let mut new_m = Array2::from_elem(
            (
                (self.motif.len() as isize + left_incr + right_incr) as usize,
                DNAMotif::MONO_CT,
            ),
            1.0 / DNAMotif::MONO_CT as f32,
        );
        let old_offset = if left_incr == -1 { 1 } else { 0 };
        let new_offset = if left_incr == 1 { 1 } else { 0 };
        let to = if right_incr == -1 {
            self.motif.len() - 1
        } else {
            self.motif.len()
        } - old_offset;
        {
            let scores = self.motif.get_scores();
            for i in 0..to {
                for b in (0..DNAMotif::MONO_CT) {
                    new_m[[new_offset + i, b]] = scores[[old_offset + i, b]];
                }
            }
        }
        self.motif = new_m.into();

        true
    }
}

pub trait MatrixPlus<'a> {
    fn eval_seqs(
        &self,
        pool: &mut Pool,
        seqs: impl Iterator<Item = &'a [u8]>,
    ) -> Vec<(&'a [u8], ScoredPos)>;
    fn normalize_scores(&mut self);
}

impl<'a> MatrixPlus<'a> for DNAMotif {
    /// apply motif to sequences in a FASTA file, returning sequences and scores
    fn eval_seqs(
        &self,
        pool: &mut Pool,
        seqs: impl Iterator<Item = &'a [u8]>,
    ) -> Vec<(&'a [u8], ScoredPos)> {
        // FIXME: b/c we wind up re-analyzing these files again and again,
        // we should probably just read into memory once and be done w/ it

        let mut results: Vec<(&'a [u8], ScoredPos)> =
            seqs.map(|seq| (seq, ScoredPos::default())).collect();
        results
            .split_iter_mut()
            .for_each(&pool.spawner(), |p| match self.score(p.0) {
                Ok(sp) => {
                    p.0 = &p.0[sp.loc..sp.loc + self.scores.dim().0];
                    p.1 = sp;
                }
                _ => (),
            });

        results
    }

    /// normalize scores in-place by summing each column and dividing each value
    fn normalize_scores(&mut self) {
        let width = self.scores.dim().0;

        for i in 0..width {
            let mut tot = 0.0;
            for j in 0..4 {
                tot += self.scores[[i, j]];
            }
            for j in 0..4 {
                self.scores[[i, j]] = self.scores[[i, j]] / tot;
            }
        }
    }
}

/// initializes array to 1.0, then increments for each base match
fn tally<'a, I, M>(motif: &M, seqs: I) -> Array2<f32>
where
    I: Iterator<Item = &'a [u8]>,
    M: Motif,
{
    let len = motif.len();
    let mut new_scores = Array2::from_elem((len, 4), 1.0);
    for ref seq in seqs {
        let loc = motif.score(*seq).expect("refine_mean - loc").loc;
        for i in 0..len {
            new_scores[[i, M::lookup(seq[loc + i]).expect("tally")]] += 1.0;
        }
    }
    new_scores
}

/// choose samples for EM
/// filter "pos" seqs by a threshold, collecting passing seqs into Vec passing_pos
/// sort "negs" seqs by score and (1) return an equal-len Vec (repr_neg) of the highest scoring "neg" seqs, and
///   find the position of the first failing seq (cutoff_neg)
pub fn choose<'a>(
    threshold: f32,
    pos_v: &mut Vec<(&'a [u8], ScoredPos)>,
    neg_v: &mut Vec<(&'a [u8], ScoredPos)>,
) -> (
    Vec<(&'a [u8], ScoredPos)>,
    usize,
    Vec<(&'a [u8], ScoredPos)>,
) {
    let passing_pos: Vec<(&'a [u8], ScoredPos)> = pos_v
        .iter()
        .filter_map(|ref t| {
            if t.1.sum >= threshold {
                Some(t.clone())
            } else {
                None
            }
        })
        .cloned()
        .collect();
    neg_v.sort_by(|&(_, ref score_a), &(_, ref score_b)| {
        score_b.sum.partial_cmp(&score_a.sum).expect("float sort")
    });
    let repr_neg = neg_v.iter().take(passing_pos.len()).cloned().collect();
    let cutoff_neg = match neg_v
        .iter()
        .enumerate()
        .find(|&(_, &(_, ref score))| score.sum <= threshold)
    {
        Some((i, _)) => i,
        None => 0,
    };

    (passing_pos, cutoff_neg, repr_neg)
}

pub fn read_seqs(fname: &str) -> Vec<Vec<u8>> {
    fasta::Reader::from_file(fname)
        .expect(format!("couldn't open {}", fname).as_str())
        .records()
        .map(|_rec| _rec.expect("unwrap record").seq().to_vec())
        .collect()
}

/// wrapper to make seqs_to_dyad more usable
pub enum SeqRefList<'a, 'b> {
    WithSP(&'b mut dyn Iterator<Item = &'a (&'a [u8], ScoredPos)>),
    Without(&'b mut dyn Iterator<Item = &'a [u8]>),
}

/// given motif, choose matching sequences.  resulting dyad has two copies of the motif (init and motif)
pub fn seqs_to_dyad<'a, 'b>(
    mut cpu_pool: &mut Pool,
    init: &DNAMotif,
    pos: SeqRefList<'a, 'b>,
    pos_ct: usize,
    neg: SeqRefList<'a, 'b>,
    neg_ct: usize,
    threshold: Option<f32>,
) -> (f64, DyadMotif<'a, DNAMotif>) {
    let mut pos_v = match pos {
        SeqRefList::WithSP(it) => init.eval_seqs(cpu_pool, it.map(|t| t.0)),
        SeqRefList::Without(it) => init.eval_seqs(cpu_pool, it),
    };
    let mut neg_v = match neg {
        SeqRefList::WithSP(it) => init.eval_seqs(cpu_pool, it.map(|t| t.0)),
        SeqRefList::Without(it) => init.eval_seqs(cpu_pool, it),
    };

    let thresh = threshold.unwrap_or_else(|| passing_threshold(init));
    let (pos_seqs, neg_hits, neg_seqs) = choose(thresh, &mut pos_v, &mut neg_v);
    (
        scaled_fisher(pos_v.len(), pos_ct, neg_hits, neg_ct),
        DyadMotif {
            init: init.clone(),
            history: vec![MotifHistory::Init],
            motif: init.clone(),
            kmer_len: (0, 0),
            gap_len: 0,
            pos_seqs: pos_seqs,
            neg_seqs: neg_seqs,
            neg_seq_ct: neg_ct,
            score: f64::NAN,
        },
    )
}

// apply motif to find sequences; take mean; return new dyad + distance calculation
pub fn em_cycle<'a>(
    mut cpu_pool: &mut Pool,
    init: &DNAMotif,
    pos: &'a Vec<(&'a [u8], ScoredPos)>,
    pos_ct: usize,
    neg: &'a Vec<(&'a [u8], ScoredPos)>,
    neg_ct: usize,
    threshold: f32,
) -> (f32, f64, DNAMotif) {
    let (p_val, dyad) = seqs_to_dyad(
        &mut cpu_pool,
        init,
        SeqRefList::WithSP(&mut pos.iter()),
        pos_ct,
        SeqRefList::WithSP(&mut neg.iter()),
        neg_ct,
        Some(threshold),
    );
    info!(
        "~~ seqs_to_dyad -> info_content={}, pos_seqs={}, threshold={}",
        dyad.motif.info_content(),
        dyad.pos_seqs.len(),
        threshold
    );
    let mean = dyad.refine_mean();

    (
        dyad.motif.distance(&mean.motif).expect("distance calc"),
        p_val,
        mean.motif,
    )
}

/// Expectation maximization (EM) algorithm which (1) finds reads matching above some threshold,
/// then (2) uses them to construct a new motif by tallying their bases.
/// Returns when another round of the EM produces the same motif.
/// Sequence counts are explcitely passed in order to support indexing, ie, passing in a subset
/// of sequences.
pub fn mean_until_stable<'a>(
    mut cpu_pool: &mut Pool,
    init: &DNAMotif,
    pos: &'a Vec<(&'a [u8], ScoredPos)>,
    pos_ct: usize,
    neg: &'a Vec<(&'a [u8], ScoredPos)>,
    neg_ct: usize,
    threshold: f32,
) -> (f64, DyadMotif<'a, DNAMotif>) {
    let (dist, p_val, motif) = em_cycle(&mut cpu_pool, init, &pos, pos_ct, &neg, neg_ct, threshold);
    info!("mean_until_stable: dist={}, EPISOLON={}", dist, EPSILON);
    if dist < EPSILON {
        info!("~~ found thing");
        let mut p = pos.iter();
        let mut n = neg.iter();
        seqs_to_dyad(
            &mut cpu_pool,
            &motif,
            SeqRefList::WithSP(&mut p),
            pos_ct,
            SeqRefList::WithSP(&mut n),
            neg_ct,
            Some(threshold),
        )
    } else {
        mean_until_stable(&mut cpu_pool, &motif, pos, pos_ct, neg, neg_ct, threshold)
    }
}

/// given a motif, calculate the passing threshold based on a simple function of information content
/// see README for more info
pub fn passing_threshold(motif: &DNAMotif) -> f32 {
    72.54231453 / (motif.info_content() + 75.48262684)
}

/// returns list of sequence indices matching motif by applying the motif to all kmers
/// in the index, and filtering by threshhold
pub fn kmer_heur(
    motif: &DNAMotif,
    kmers: &KmerIndex,
    exclude: Option<&Range<usize>>,
) -> Result<HashSet<usize>, PSSMError> {
    let threshold: f32 = passing_threshold(motif);

    let mut all_ids = hashset![];
    for (kmer, ids) in &kmers.0 {
        match motif.score(&kmer[..]) {
            Err(PSSMError::NullMotif) => {
                // ignore degenerate edge-cases
            }
            Err(e) => return Err(e),
            Ok(sp) => {
                if sp.sum >= threshold {
                    if let Some(r) = exclude {
                        // ignore matches overlapping with excluded range
                        all_ids.extend(ids.iter().filter_map(|i| {
                            if r.contains(i) {
                                None
                            } else {
                                Some(*i)
                            }
                        }));
                    } else {
                        all_ids.extend(ids.iter());
                    }
                }
            }
        }
    }
    Ok(all_ids)
}

#[cfg(test)]
mod tests {
    use super::*;
    const MOTIF: &'static [u8] = b"GGCCTAGCCATG";
    const POS_FNAME: &'static str = "pos.fa"; // all contain MOTIF
                                              //const POS_FNAME: &'static str = "test.fa"; // various motifs at various frequencies
    const NEG_FNAME: &'static str = "neg.fa";

    #[test]
    #[ignore]
    fn kmers_to_m() {
        let m = DyadMotif::<DNAMotif>::kmers_to_matrix(b"ATGC", 1, b"ATGC");
        let expected = Array::from_vec(vec![
            0.85, 0.05, 0.05, 0.05, 0.05, 0.85, 0.05, 0.05, 0.05, 0.05, 0.85, 0.05, 0.05, 0.05,
            0.05, 0.85, 0.25, 0.25, 0.25, 0.25, 0.85, 0.05, 0.05, 0.05, 0.05, 0.85, 0.05, 0.05,
            0.05, 0.05, 0.85, 0.05, 0.05, 0.05, 0.05, 0.85,
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

    fn choose_sm(
        pos_v: &mut Vec<(Vec<u8>, f64)>,
        neg_v: &mut Vec<(Vec<u8>, f64)>,
    ) -> Option<(Vec<Vec<u8>>, Vec<Vec<u8>>)> {
        pos_v.sort_by(|&(_, score_a), &(_, score_b)| {
            score_b.partial_cmp(&score_a).expect("float sort")
        });
        neg_v.sort_by(|&(_, score_a), &(_, score_b)| {
            score_b.partial_cmp(&score_a).expect("float sort")
        });
        Some((
            pos_v
                .iter()
                .map(|&(ref s, _)| s.clone())
                .take(100)
                .collect(),
            neg_v
                .iter()
                .map(|&(ref s, _)| s.clone())
                .take(100)
                .collect(),
        ))
    }

    #[test]
    #[ignore]
    fn test_find_one_motif() {
        println!("dyad::test_find");
        let v = DyadMotif::<DNAMotif>::passing_kmers(POS_FNAME, NEG_FNAME, vec![].as_ref());
        let pos = read_seqs(POS_FNAME);
        let neg = read_seqs(NEG_FNAME);
        let dyads = DyadMotif::<DNAMotif>::motifs(v, &pos, &neg, choose);
        //dyads[0].refine(100);
    }

    #[test]
    fn test_info_content() {
        let m = DNAMotif::from(
            Array::from_vec(vec![
                1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            ])
            .into_shape((3, 4))
            .unwrap(),
        );
        assert_eq!(m.info_content(), 6.0);

        let m = DNAMotif::from(
            Array::from_vec(vec![
                0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
            ])
            .into_shape((3, 4))
            .unwrap(),
        );
        assert_eq!(m.info_content(), 0.0);
    }

    #[test]
    fn test_find_motifs() {
        env_logger::init();
        info!("################### test->test_find_motifs ###################");
        let v = DyadMotif::<DNAMotif>::passing_kmers(POS_FNAME, NEG_FNAME, vec![].as_ref());
        info!("{} passing kmers", v.len());
        //find_motifs(v, POS_FNAME, NEG_FNAME);
    }

    #[test]
    #[ignore]
    fn blah() {
        env_logger::init();

        info!("################### test->blah ###################");

        let m: Array2<f32> = Array::from(vec![
            0.9962687,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.9962687,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.9962687,
            0.0012437812,
            0.9962687,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.9962687,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.9962687,
            0.0012437812,
            0.9962687,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.0012437811,
            0.0012437811,
            0.0012437811,
            0.99626863,
            0.9962687,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.9962687,
            0.0012437812,
            0.0012437812,
            0.0012437812,
            0.9962687,
            0.0012437812,
            0.0012437812,
        ])
        .into_shape((11, 4))
        .unwrap();
        let motif = DNAMotif::from(m);
        println!("score: {:?}", motif.score(b"AGGGCAAGTAGCTGATTGAAGTAGCAGAGGCTGGCGTCCAAGCGGTAATAAACAAGCGATGAAAAATACTGAAGTACCTGGAGCCATTACTCAATAGGAGCAGTCCGTGAAACCTGTGCGGCGTGTAGCGAATGTTCGGCACATTATGTAAGAAGACATTTGCTTATTCACGAAATCAGCGCGAACCCTCCATCTGCCGA".iter()));
        let mut pool = make_pool(*CPU_COUNT).unwrap();

        let pos_seqs = read_seqs(POS_FNAME);
        let neg_seqs = read_seqs(NEG_FNAME);

        let mut p = motif.eval_seqs(&mut pool, pos_seqs.iter().map(|s| s.as_slice()));
        let mut n = motif.eval_seqs(&mut pool, neg_seqs.iter().map(|s| s.as_slice()));
        println!("p: {:?}", p.iter().map(|t| t.1.sum).collect::<Vec<f32>>());
        let (pos_seqs, _, neg_seqs) = choose(0.97, &mut p, &mut n);
        println!("pos_seqs.len: {}", pos_seqs.len());
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
    fn test_slide() {
        let m1 = DNAMotif::from_degenerate(b"NAAAA").unwrap();
        let mut d1 = DyadMotif {
            init: m1.clone(),
            history: vec![MotifHistory::Init],
            motif: m1.clone(),
            kmer_len: (0,0),
            gap_len: 0,
            pos_seqs: vec![],
            neg_seqs: vec![],
            neg_seq_ct: 0,
            score: 0.0,
        };
        assert!(d1.slide());
        println!(
            "-- new degen: {}",
            str::from_utf8(d1.motif.degenerate_consensus().as_slice()).unwrap()
        );
        assert_eq!(d1.motif.degenerate_consensus(), b"AAAAN".to_vec());
        assert_eq!(
            d1.history,
            vec![
                MotifHistory::Init,
                MotifHistory::ContractLeft,
                MotifHistory::ExpandRight
            ]
        );

        let m2 = DNAMotif::from_degenerate(b"AAAAN").unwrap();
        let mut d2 = DyadMotif {
            init: m2.clone(),
            history: vec![MotifHistory::Init],
            motif: m2.clone(),
            kmer_len: (0,0),
            gap_len: 0,
            pos_seqs: vec![],
            neg_seqs: vec![],
            neg_seq_ct: 0,
            score: 0.0,
        };

        assert!(d2.slide());
        assert_eq!(d2.motif.degenerate_consensus(), b"NAAAA".to_vec());
        assert_eq!(
            d2.history,
            vec![
                MotifHistory::Init,
                MotifHistory::ExpandLeft,
                MotifHistory::ContractRight
            ]
        );
    }
}
