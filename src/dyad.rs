use darwin_rs::{Individual, Population, PopulationBuilder, SimulationBuilder};
use jobsteal::{make_pool, BorrowSpliterator, IntoSpliterator, Pool, Spliterator};
use std::cmp::{max, min};
use std::str;
use std::{f64, usize};
//use darwin_rs::select::MaximizeSelector;
use fishers_exact::{fishers_exact, TestTails};

use bio::io::fasta;
use bio::pattern_matching::pssm::{Motif, ScoredPos};
use ndarray::prelude::{Array, Array2};
use rand;
use rand::Rng;

use super::*;
use ctr::*;

const P_CUTOFF: f64 = 0.001;
const MODAL_MAX_SEQS: usize = 400;

#[derive(Debug, Clone)]
pub enum MotifHistory {
    Init,
    Mutate,
    Mean,
    MeanDiff,
    Mode,
    Reset,
}

#[derive(Debug, Clone)]
pub struct DyadMotif<M>
where
    M: Motif + Clone,
{
    /// initial state based on kmers
    init: M,
    /// history of actions that produced this motif
    pub history: Vec<MotifHistory>,
    /// weights updated by GA
    pub motif: M,
    /// kmer len
    pub kmer_len: usize,
    /// gap len
    pub gap_len: usize,
    /// sequences matching our motif
    pub pos_seqs: Vec<(Vec<u8>, ScoredPos)>,
    /// sequences representing background
    pub neg_seqs: Vec<(Vec<u8>, ScoredPos)>,
    /// score - sum(neg) / sum(pos)
    pub score: f64,
}

impl<M> DyadMotif<M>
where
    M: Motif + Clone + Sync + Send + From<Array2<f32>>,
{
    fn fasta_to_ctr(fname: &str) -> (GappedKmerCtr<M>, usize) {
        let mut ctr = GappedKmerCtr::new(KMER_LEN, MIN_GAP, MAX_GAP);
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

    fn hamming(a: &[u8], b: &[u8]) -> usize {
        a.iter()
            .zip(b.iter())
            .map(|(a, b)| if *a == *b { 1 } else { 0 })
            .sum()
    }

    /// P-values returned by the Fisher exact test don't change meaningfully as the values get larger.
    /// eg, both [100, 200, 10_000, 10_000] and [1_000, 2_000, 100_000, 100_000] yield a P-value well
    /// below our cutoff.  therefore, we can safely scale the values down if they're above some arbitrary
    /// threshold.
    fn scaled_fisher(_ct1: usize, _tot1: usize, _ct2: usize, _tot2: usize) -> f64 {
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

    /// generate kmers, tablulate, and apply Fisher exact test
    pub fn passing_kmers(pos_fname: &str, neg_fname: &str) -> Vec<(usize, usize, usize, f64)> {
        let (pos, pos_ct) = DyadMotif::<DNAMotif>::fasta_to_ctr(pos_fname);
        let (neg, neg_ct) = DyadMotif::<DNAMotif>::fasta_to_ctr(neg_fname);

        let (width, height, gap) = pos.ctr.dim();

        let mut pool = make_pool(*CPU_COUNT).unwrap();
        info!("passing_kmers - created pool of {} threads", *CPU_COUNT);
        let mut indices: Vec<(usize, Vec<(usize, usize, f64)>)> =
            (0..width).map(|i| (i, vec![])).collect();
        indices
            .split_iter_mut()
            .for_each(&pool.spawner(), |&mut (ref i, ref mut v)| {
                //for i in 0..width {
                for j in 0..height {
                    for k in 0..gap {
                        if pos.ctr[[*i, j, k]] > neg.ctr[[*i, j, k]] {
                            let p = DyadMotif::<DNAMotif>::scaled_fisher(
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
        let mut passing = Vec::new();
        for (i, v) in indices {
            passing.extend(v.iter().map(|&(j, k, p)| (i, j, k, p)));
        }
        passing
    }

    pub fn motifs<F>(
        chosen: Vec<(usize, usize, usize, f64)>,
        pos_fname: &str,
        neg_fname: &str,
        chooser: F,
    ) -> Vec<DyadMotif<DNAMotif>>
    where
        F: Fn(
            &mut Vec<(Vec<u8>, ScoredPos)>,
            &mut Vec<(Vec<u8>, ScoredPos)>,
        ) -> Option<(Vec<(Vec<u8>, ScoredPos)>, Vec<(Vec<u8>, ScoredPos)>)>,
    {
        info!("using {} cpus", *CPU_COUNT);
        let mut pool = make_pool(*CPU_COUNT).unwrap();
        let mut dyads = Vec::new();
        for (idx, &(i, j, k, _)) in chosen.iter().enumerate() {
            if idx % 500 == 0 {
                info!("creating dyad #{} / {}", idx, chosen.len());
            }

            let init: DNAMotif = (DyadMotif::<DNAMotif>::kmers_to_matrix(
                GappedKmerCtr::<DNAMotif>::int_to_kmer(KMER_LEN, i).as_slice(),
                k,
                GappedKmerCtr::<DNAMotif>::int_to_kmer(KMER_LEN, j).as_slice(),
            ))
            .into();

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

            let mut pos_v = init.eval_file(&mut pool, pos_fname);
            let mut neg_v = init.eval_file(&mut pool, neg_fname);
            let (pos_seqs, neg_seqs) =
                chooser(&mut pos_v, &mut neg_v).expect("motifs found bad one (1)");

            info!(
                "DyadMotif::motifs - {} / {} seqs used",
                pos_seqs.len(),
                pos_v.len()
            );

            let copy = init.clone();

            let mut d = DyadMotif {
                init: init,
                history: vec![MotifHistory::Init],
                motif: copy,
                kmer_len: KMER_LEN,
                gap_len: MIN_GAP + k,
                pos_seqs: pos_seqs,
                neg_seqs: neg_seqs,
                score: f64::NAN,
            };
            d.calculate_fitness();
            dyads.push(d);
        }
        dyads
    }

    /// use a genetic algorithm - generate @mut_ct mutants, mutate further, and return the best
    pub fn refine_GA(&self, mut_ct: usize) -> DyadMotif<M> {
        // make an initial population of 100 copies of the motif
        let mut init_pop = (0..mut_ct)
            .map(|_| self.clone())
            .collect::<Vec<DyadMotif<M>>>();
        for ind in init_pop.iter_mut() {
            ind.mutate();
            ind.history.push(MotifHistory::Mutate);
        }

        let means_based = self.clone().refine_mean();
        init_pop.push(means_based);

        let meansdiv_based = self.clone().refine_meandiv();
        init_pop.push(meansdiv_based);

        let mode_based = self.clone().refine_via_mode();
        init_pop.push(mode_based);

        let population1 = PopulationBuilder::<DyadMotif<M>>::new()
            .initial_population(init_pop.as_slice())
            .increasing_exp_mutation_rate(1.03)
            .reset_limit_end(0) // disable resetting
            .finalize()
            .expect("PopulationBuilder");
        let mut sim = SimulationBuilder::<DyadMotif<M>>::new()
            .iterations(11)
            .threads(1)
            .add_population(population1)
            .finalize()
            .expect("some problem making builder");

        sim.run();
        sim.print_fitness();
        //sim.simulation_result.fittest[0].individual.calculate_fitness();

        sim.simulation_result.fittest[0].individual.clone()
    }

    /// simple means-based refinement
    pub fn refine_mean(&self) -> DyadMotif<M> {
        let len = self.motif.len();
        let new_scores = tally(&self.motif, self.pos_seqs.iter().map(|&(ref seq, _)| seq));

        let m: M = new_scores.into();
        //m.normalize_scores();
        //m.calc_minmax();

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
            score: f64::NAN,
        };

        d.calculate_fitness();
        d
    }

    /// tally base-matches for positive seqs as in refine_mean, but also tally
    /// negetive seqs and divide positive tallies by negative before normalizing
    pub fn refine_meandiv(&self) -> DyadMotif<M> {
        let pos_tally = tally(&self.motif, self.pos_seqs.iter().map(|&(ref seq, _)| seq));
        let neg_tally = tally(&self.motif, self.neg_seqs.iter().map(|&(ref seq, _)| seq));

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
            score: f64::NAN,
        };

        d.calculate_fitness();
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

                diffs[[i, j]] = Self::hamming(a, b);
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
        //motif.normalize_scores();
        //motif.calc_minmax();
        info!(
            "mode: {:?}",
            String::from_utf8(motif.degenerate_consensus()).unwrap()
        );

        let mut pos: Vec<(Vec<u8>, ScoredPos)> = self
            .pos_seqs
            .iter()
            .map(|&(ref seq, _)| (seq.clone(), motif.score(seq).expect("mode pos")))
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

        let mut neg: Vec<(Vec<u8>, ScoredPos)> = self
            .neg_seqs
            .iter()
            .map(|&(ref seq, _)| (seq.clone(), motif.score(seq).expect("mode neg")))
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
            score: f64::NAN,
        };
        d.calculate_fitness();
        d
    }

    pub fn refine(&self, mut_ct: usize) -> DyadMotif<M> {
        self.refine_GA(mut_ct)
    }

    /// stringified self.motif.degenerate_consensus()
    pub fn show_motif(&self) -> String {
        String::from_utf8(self.motif.degenerate_consensus()).expect("show_motif")
    }
}

pub trait MatrixPlus {
    fn eval_file(&self, pool: &mut Pool, fname: &str) -> Vec<(Vec<u8>, ScoredPos)>;
    fn normalize_scores(&mut self);
}

impl MatrixPlus for DNAMotif {
    /// apply motif to sequences in a FASTA file, returning sequences and scores
    fn eval_file(&self, pool: &mut Pool, fname: &str) -> Vec<(Vec<u8>, ScoredPos)> {
        // FIXME: b/c we wind up re-analyzing these files again and again,
        // we should probably just read into memory once and be done w/ it
        let mut v = Vec::new();
        for _rec in fasta::Reader::from_file(fname)
            .expect(format!("couldn't open {}", fname).as_str())
            .records()
        {
            let rec = _rec.expect("unwrap record");
            v.push((rec.seq().to_vec(), ScoredPos::default()));
        }

        if v.len() == 0 {
            panic!("empty file: {}", fname);
        }

        v.split_iter_mut()
            .for_each(&pool.spawner(), |p| match self.score(&p.0) {
                Ok(sp) => {
                    p.1 = sp;
                }
                _ => (),
            });
        v
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
    I: Iterator<Item = &'a Vec<u8>>,
    M: Motif,
{
    let len = motif.len();
    let mut new_scores = Array2::from_elem((len, 4), 1.0);
    for ref seq in seqs {
        let loc = motif.score(seq.as_slice()).expect("refine_mean - loc").loc;
        for i in 0..len {
            new_scores[[i, M::lookup(seq[loc + i]).expect("tally")]] += 1.0;
        }
    }
    new_scores
}

/// choose samples for EM
pub fn choose(
    pos_v: &mut Vec<(Vec<u8>, ScoredPos)>,
    neg_v: &mut Vec<(Vec<u8>, ScoredPos)>,
) -> Option<(Vec<(Vec<u8>, ScoredPos)>, Vec<(Vec<u8>, ScoredPos)>)> {
    pos_v.sort_by(|&(_, ref score_a), &(_, ref score_b)| {
        score_b.sum.partial_cmp(&score_a.sum).expect("float sort")
    });
    neg_v.sort_by(|&(_, ref score_a), &(_, ref score_b)| {
        score_b.sum.partial_cmp(&score_a.sum).expect("float sort")
    });

    let mut cutoff = 0;
    for (i, &(_, ref score)) in pos_v.iter().enumerate() {
        if score.sum <= 0.9 {
            break;
        }
        cutoff = i;
    }
    if cutoff == 0 {
        warn!(
            "-- bad: pos={:?}/{}, {:?}/{}, neg={:?}/{}, {:?}/{}",
            str::from_utf8(&pos_v[0].0).unwrap(),
            pos_v[0].1.sum,
            str::from_utf8(&pos_v[1].0).unwrap(),
            pos_v[1].1.sum,
            str::from_utf8(&neg_v[0].0).unwrap(),
            neg_v[0].1.sum,
            str::from_utf8(&neg_v[1].0).unwrap(),
            neg_v[1].1.sum,
        );
        return None;
    }
    Some((
        pos_v.iter().take(cutoff).cloned().collect(),
        neg_v.iter().take(cutoff).cloned().collect(),
    ))
}
/*
fn crossover_motifs(
    mine: &mut Motif,
    theirs: &mut Motif,
    pos_fname: &str,
    neg_fname: &str,
) -> Motif {
    assert_eq!(mine.len(), theirs.len());

    // store sequences in memory, as IO doesn't play nice w/ parallelism
    let mut pos_seqs = Vec::new();
    for _rec in fasta::Reader::from_file(pos_fname)
        .expect(format!("couldn't open {}", pos_fname).as_str())
        .records()
    {
        let rec = _rec.expect("unwrap record");
        pos_seqs.push((rec.seq().to_vec(), ScoredPos::default(), ScoredPos::default()));
    }
    let mut neg_seqs = Vec::new();
    for _rec in fasta::Reader::from_file(neg_fname)
        .expect(format!("couldn't open {}", neg_fname).as_str())
        .records()
    {
        let rec = _rec.expect("unwrap record");
        neg_seqs.push((rec.seq().to_vec(), ScoredPos::nil(), ScoredPos::nil()));
    }

    // step one: reduce input to match the width of motif
    // this necessarily means chosing between the best match position from @mine or @theirs
    fn max_slice(
        pwm_len: usize,
        mine: &Motif,
        theirs: &Motif,
        seq: &[u8],
    ) -> (ScoredPos, ScoredPos) {
        let s = mine.score(seq).expect(
            format!("self couldn't score: {:?}", seq)
                .as_str(),
        );
        let o = theirs.score(seq).expect(
            format!("other couldn't score: {:?}", seq)
                .as_str(),
        );

        let t = if s.loc == o.loc {
            (s, o)
        } else if o.sum > s.sum {
            let x = mine.score(&seq[o.loc..o.loc + pwm_len]).expect(
                "couldn't score slice (1)",
            );
            (x, o)
        } else {
            let x = theirs.score(&seq[s.loc..s.loc + pwm_len]).expect(
                "couldn't score slice (2)",
            );
            (s, x)
        };

        t
    }

    let mut pool = make_pool(*CPU_COUNT).unwrap();
    let pwm_len = mine.len();
    pos_seqs.split_iter_mut().for_each(&pool.spawner(), |t| {
        let scores = max_slice(pwm_len, mine, theirs, t.0.as_slice());
        t.1 = scores.0;
        t.2 = scores.1;
    });
    neg_seqs.split_iter_mut().for_each(&pool.spawner(), |t| {
        let scores = max_slice(pwm_len, mine, theirs, t.0.as_slice());
        t.1 = scores.0;
        t.2 = scores.1;
    });

    // step 2: create new PWM by choosing the best base at each position
    let mut new_m = Array2::zeros((pwm_len, 4));
    for i in 0..pwm_len {
        // positive set
        let mut s_ptally = 0.0;
        let mut o_ptally = 0.0;
        for &(_, ref sscore, ref oscore) in pos_seqs.iter() {
            s_ptally += sscore.scores[i];
            o_ptally += oscore.scores[i];
        }

        // negative set
        let mut s_ntally = 0.0;
        let mut o_ntally = 0.0;
        for &(_, ref sscore, ref oscore) in neg_seqs.iter() {
            s_ntally += sscore.scores[i];
            o_ntally += oscore.scores[i];
        }
        //println!("@ i={}, o > s? {:?}", i, (o_ptally / o_ntally) > (s_ptally / s_ntally));
        let mut __c: usize = 0;
        for b in 0..4 {
            new_m[[i, b]] = if (o_ptally / o_ntally) > (s_ptally / s_ntally) {
                __c += 1;
                theirs.scores[[i, b]]
            } else {
                mine.scores[[i, b]]
            }
        }
        if __c != 0 && __c != 4 {
            println!("@@ weird mixed");
        }
    }

    Motif::from(new_m)
}
*/

impl<M> Individual for DyadMotif<M>
where
    M: Motif + Clone + Sync + Send + From<Array2<f32>>,
{
    //const CAN_CROSSOVER: bool = false;

    /// shift value at each position
    fn mutate(&mut self) {
        // Mutate the scores
        let mut scores = self.motif.get_scores().clone();
        for x in scores.iter_mut() {
            // r is on (0,1)
            let r = rand::random::<f32>();
            // by subtracting 0.5, we allow for negative random numbers
            // by scaling by 0.02, we limit changes to (-0.01,0.01)
            let new_x = *x + MUT_INCR * (r - 0.5);
            *x = if new_x < 0.0 { 0.0 } else { new_x };
        }
        self.motif = scores.into();
    }

    fn calculate_fitness(&mut self) -> f64 {
        // Calculate how good the data values are compared to the perfect solution

        let mut pool = make_pool(*CPU_COUNT).unwrap();

        let pos_sum: f64 = self
            .pos_seqs
            .clone()
            .split_iter()
            .map(|&(ref seq, _)| self.motif.score(seq.iter()).expect("score?").sum as f64)
            .collect::<Vec<f64>>(&pool.spawner())
            .iter()
            .sum();

        let neg_sum: f64 = self
            .neg_seqs
            .clone()
            .split_iter()
            .map(|&(ref seq, _)| self.motif.score(seq.iter()).expect("score?").sum as f64)
            .collect::<Vec<f64>>(&pool.spawner())
            .iter()
            .sum();

        let score = if pos_sum == 0.0 {
            f64::INFINITY
        } else if neg_sum == 0.0 {
            0.0
        } else {
            if self.pos_seqs.len() == self.neg_seqs.len() {
                neg_sum / pos_sum
            } else {
                let pos: f64 = pos_sum / self.pos_seqs.len() as f64;
                let neg: f64 = neg_sum / self.neg_seqs.len() as f64;

                neg / pos
            }
        };
        self.score = score;
        score
    }

    /// initialize array with random values, then normalize
    /// so each position sums to 1.0
    fn reset(&mut self) {
        println!("-- reset");
        // bases == 4
        self.motif = self.init.clone();
        self.history.push(MotifHistory::Reset);
        self.score = f64::NAN;
    }

    /* FIXME: switched from crossover fork to regular darwin-rs

    fn crossover(&mut self, other: &mut Self) -> Self {
        info!("DyadMotif::crossover");
        let new_motif = crossover_motifs(
            &mut self.motif,
            &mut other.motif,
            self.pos_fname,
            self.neg_fname,
        );

        DyadMotif {
            init: self.motif.clone(),
            motif: new_motif,
            ctr: self.ctr.clone(),
            pos_fname: self.pos_fname,
            neg_fname: self.neg_fname,
        }
        self.clone()
    }*/
}

pub fn find_motifs(
    chosen: Vec<(usize, usize, usize, f64)>,
    pos_fname: &str,
    neg_fname: &str,
) -> Vec<DyadMotif<DNAMotif>> {
    let mut pool = make_pool(*CPU_COUNT).unwrap();
    let motifs = DyadMotif::<DNAMotif>::motifs(chosen, pos_fname, neg_fname, choose);
    motifs
        .iter()
        .enumerate()
        .map(|(_, dyad)| {
            // dyad wraps the sequences chosen by our method
            let mut new_dyad = dyad.refine(100);

            // now that the GA has [hopefully] improved our PWM, we need to choose new seqs
            let mut pos_v = new_dyad.motif.eval_file(&mut pool, pos_fname);
            let mut neg_v = new_dyad.motif.eval_file(&mut pool, neg_fname);
            let (pos_seqs, neg_seqs) =
                choose(&mut pos_v, &mut neg_v).expect("motifs found bad one (2)");
            new_dyad.pos_seqs = pos_seqs;
            new_dyad.neg_seqs = neg_seqs;

            //new_dyad.repopulate_seqs(pos_fname, neg_fname, choose);
            let new2 = new_dyad.refine(100);
            new2
        })
        .collect()
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
        let v = DyadMotif::<DNAMotif>::passing_kmers(POS_FNAME, NEG_FNAME);
        let dyads = DyadMotif::<DNAMotif>::motifs(v, POS_FNAME, NEG_FNAME, choose);
        dyads[0].refine(100);
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
        let v = DyadMotif::<DNAMotif>::passing_kmers(POS_FNAME, NEG_FNAME);
        info!("{} passing kmers", v.len());
        find_motifs(v, POS_FNAME, NEG_FNAME);
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

        let mut p = motif.eval_file(&mut pool, POS_FNAME);
        let mut n = motif.eval_file(&mut pool, NEG_FNAME);
        println!("p: {:?}", p.iter().map(|t| t.1.sum).collect::<Vec<f32>>());
        let (pos_seqs, neg_seqs) = choose(&mut p, &mut n).expect("motifs found bad one (2)");
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
}
