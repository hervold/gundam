use super::*;
use jobsteal::{make_pool, BorrowSpliterator, IntoSpliterator, Pool, Spliterator};
use std::collections::HashMap;
use std::collections::HashSet;

const KMER_LEN: usize = 8;
const KMER_INCR: usize = 4;

pub type Kmer = [u8; KMER_LEN];

#[derive(Debug)]
pub struct KmerIndex(pub HashMap<Kmer, HashSet<usize>>);

impl KmerIndex {
    pub fn new(seqs: &Vec<Vec<u8>>) -> KmerIndex {
        let mut kmers_to_idx = HashMap::new();
        for (seq_i, seq) in seqs.iter().enumerate() {
            let mut start = 0;
            while start + KMER_LEN < seq.len() {
                let mut kmer = [0u8; KMER_LEN];
                for i in 0..KMER_LEN {
                    kmer[i] = seq[start + i];
                }
                start += KMER_INCR;

                kmers_to_idx
                    .entry(kmer)
                    .or_insert_with(HashSet::new)
                    .insert(seq_i);
            }
        }
        KmerIndex(kmers_to_idx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ergo_std::*;

    #[test]
    fn test_constructor() {
        let kmer1: Kmer = [b'A', b'T', b'G', b'C', b'A', b'T'];
        let kmer2: Kmer = [b'T', b'G', b'C', b'A', b'T', b'G'];
        let kmer3: Kmer = [b'G', b'C', b'A', b'T', b'G', b'C'];
        let kmer4: Kmer = [b'C', b'A', b'T', b'G', b'C', b'A'];

        let kmer10: Kmer = [b'G', b'G', b'G', b'A', b'T', b'G'];
        let kmer11: Kmer = [b'A', b'T', b'G', b'C', b'G', b'G'];
        let kmer12: Kmer = [b'C', b'G', b'G', b'G', b'G', b'G'];
        let kmer13: Kmer = [b'G', b'G', b'G', b'G', b'G', b'G'];

        let seqs = vec![
            b"ATGCATGCATGCATGCATGCATGC".to_vec(),
            b"GGGATGCGGGGGGGGGGGGGGGGG".to_vec(),
        ];
        let idx = KmerIndex::new(&seqs);
        assert_eq!(
            idx.0.keys().cloned().collect::<HashSet<Kmer>>(),
            hashset![kmer1, kmer2, kmer3, kmer4, kmer10, kmer11, kmer12, kmer13]
        );

        assert_eq!(idx.0.get(&kmer1), Some(&hashset![0]));
        assert_eq!(idx.0.get(&kmer10), Some(&hashset![1]));        
        assert_eq!(idx.0.get(&[b'A'; KMER_LEN]), None);
    }
}
