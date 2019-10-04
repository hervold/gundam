# Gundam: a motif discovery algorithm based on GADEM

## Introduction

Gudam is a slapdash re-implementation of the [GADEM](https://www.ncbi.nlm.nih.gov/pubmed/19193149) algorithm, which is best described by its authors:

> Candidate words (four to six nucleotides) for constructing spaced dyads are prioritized by their degree of overrepresentation in the input sequence data. Spaced dyads are converted into starting position weight matrices (PWMs). GADEM then employs a genetic algorithm (GA), with an embedded EM algorithm to improve starting PWMs, to guide the evolution of a population of spaced dyads toward one whose entropy scores are more statistically significant

However, while GADEM was designed to find motifs enriched in a ChIP dataset, Gundam is intended to find motifs enriched in one set of sequences as compared to another.

## Overview

### K-mer selection

Gundam begins by encoding each pair of fixed-length K-mers as a pair of integers.  These serve, along with gap length, as the index of a matrix with dimensions (k, k, i), where k = 4<sup>KMER_LEN</sup> and i = the maximum gap length.  Two such matrices are computed, one from the "positive" dataset and the other from the "negative," and only dyads passing the Fisher's exact test for significance as used.  Selected kmer pairs are encoded as a position-specific scoring matrix (using the [rust-pssm library](https://github.com/ortho-the-d-stice/rust-pssm)), with matching nucleotides initially set to 0.8 and non-matching to 0.05.  Gap bases are naively assigned a weight of 0.25.  For example, given two 4-base K-mers with sequence "ATGC," and a 3-base gap between them, the following PSSM would be created:


| init:| A  | T  | G  | C  | N  | N  | N  | A  | T  | G  | C  |
| ---- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  A   |0.8 |0.05|0.05|0.05|0.25|0.25|0.25|0.8 |0.05|0.05|0.05|
|  T   |0.05|0.8 |0.05|0.05|0.25|0.25|0.25|0.05|0.8 |0.05|0.05|
|  G   |0.05|0.05|0.8 |0.05|0.25|0.25|0.25|0.05|0.05|0.8 |0.05|
|  C   |0.05|0.05|0.05|0.8 |0.25|0.25|0.25|0.05|0.05|0.05|0.8 |

### Genetic algorithm (GA) and Expectation Maximization (EM)

Gundam, much like GADEM, uses an EM algorithm to refine the PSSMs initialized by K-mer pairs.  After the PSSM is created, samples are chosen from the full set of positive and negative sequences, and a GA is used to improve PSSM weights.  The [darwin-rs](https://github.com/willi-kappler/darwin-rs) library is used for GA refinement, with fitness for each PSSM calculated as:
```
fitness = sum(negative_scores) / sum(positive_scores)
```
Note that `darwin-rs` attempts to minimize fitness scores, while PSSM matches are scored 0-1, with a perfect match scoring 1.

After a round of refinement, the best PSSM is used to score all "positive" and "negative" sequences, and all "positive" reads scoring >= 0.9 used for the next round of refinement, along with an equal number of top-scoring "negative" sequences.  This cycle is currently repeated 3 times to yield the mature PSSM.

### Usage

Currently, Gundam is divided into two separate executables in order to work around a memory leak in the [Fisher's exact test library](https://github.com/cpearce/fishers_exact) (which is really a Rust wrapper around C code extracted from the R interpreter).  In order to generate a list of K-mer pairs, use the `gundam_kmer` command:
```
gundam_kmer pos.fa neg.fa > kmer_indices
```
This should generate a file full of comma-separated integers, which is then passed into the primary executable, `gundam`:
```
gundam kmer_indices pos.fa neg.fa
```

### TODO
So many things ...

- output actual PSSMs instead of degenerate sequence representation
- longer gaps will likely require crossing-over in the GA, which isn't supported by `darwin-rs`.  I have a fix in the works.
- replace fixed cutoff score with linear model: cutoff = -0.0416887 * information_content + 0.990597
  (see Figure_1-rnd_seq_score.png - heat map with mean score regression as red line, and 2*std-dev as yellow line)



===

new scheme:
 - "degenerate" representation is canonical
 - 10-fold cross-validation:
    split sequences into 10
    motifs are generated from 9/10 and _evaluated_ on the remaining 1/10
    toss motif / eval-pval into data structure: degen -> [pvals of test sets]
    final judgement on average of pvals - constrain on freq of occurance?  need to see it 10x?

===

# PIPELINE

./target/release/gundam_kmer train_pos.fa train_neg.fa > 2019sept19.1
>> 27m50.603s

time ./target/release/gundam 2019sept19.1 train_pos.fa train_neg.fa  > 2019sept19.2

==

RUST_LOG=info ./target/release/gundam_kmer_to_degen fake_kmers.txt fake_pos.fa fake_neg.fa > fake_degen.txt
1h20m


==
no kmers; too many allocs (input size: 3) - 971.58user 304.70system 2:43.92elapsed 778%CPU (0avgtext+0avgdata 83452maxresident)
reduced allocs: 966.95user 285.15system 2:46.21elapsed 753%CPU (0avgtext+0avgdata 82248maxresident
