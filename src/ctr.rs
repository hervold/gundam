use bio::pattern_matching::pssm::{DNAMotif, Motif};
use ndarray::Array3;
use std::cmp::{max, min};
use std::convert::AsRef;
use std::marker::PhantomData;
use std::ops::Index;
use std::{mem, str, usize};

pub const MIN_KMER_LEN: usize = 4;
/// this is needed to limit memory consumption
/// FIXME: a hashmap would be better if this stuff is sparse
pub const TOT_KMER_LEN: usize = 15;

#[derive(Debug, Clone)]
pub struct GappedKmerCtr<M>
where
    M: Motif,
{
    pub ctr: Array3<usize>,
    pub kmer_len: (usize, usize),
    pub min_gap: usize,
    pub max_gap: usize,
    pub phantom: PhantomData<M>,
}

impl<M> GappedKmerCtr<M>
where
    M: Motif,
{
    /// simple bitshifting conversion of sequences to integer
    /// note that this only makes sense for fixed-width kmers, as, eg,
    /// kmer_to_int(b"A") == kmer_to_int(b"AA") == kmer_to_int(b"AAA") == 0
    pub fn kmer_to_int(kmer: &[u8]) -> usize {
        let bits_per = M::get_bits().ceil() as usize;
        if bits_per * kmer.len() > 8 * mem::size_of::<usize>() {
            panic!("kmer too large");
        }

        let mut val: usize = 0;
        for b in kmer.iter() {
            val <<= bits_per;
            val |= M::lookup(*b).expect("kmer_to_int on unknown base");
        }
        val
    }

    ///
    pub fn int_to_kmer(kmer_len: usize, i: usize) -> Vec<u8> {
        let mut kmer: Vec<u8> = Vec::new();
        let mut val = i;
        for _ in 0..kmer_len {
            kmer.push(M::rev_lk(val & 0b11));
            val >>= 2;
        }
        kmer.reverse();
        kmer
    }

    /// initalize 3d matrix for pair of equal-length kmers with a gap range, [min_gap, max_gap)
    pub fn new(kmer_len: usize, min_gap: usize, max_gap: usize) -> GappedKmerCtr<M> {
        let len: usize = 4_usize.pow(kmer_len as u32) + 1;
        GappedKmerCtr {
            ctr: Array3::zeros((len, len, max_gap - min_gap + 1)),
            kmer_len: (kmer_len, kmer_len),
            min_gap: min_gap,
            max_gap: max_gap,
            phantom: PhantomData,
        }
    }

    /// initialize a "rectangular" matrix, ie, one where kmer lengths are not equal; gap length is fixed
    pub fn rect(width: usize, gap: usize, height: usize) -> GappedKmerCtr<M> {
        let left: usize = 4_usize.pow(width as u32);
        let right: usize = 4_usize.pow(height as u32);

        debug!(
            "~~~~~~~~~~~~~~~~~~~~ allocating {} x {} x {}",
            left, right, 1
        );
        GappedKmerCtr {
            ctr: Array3::zeros((left, right, 1)),
            kmer_len: (width, height),
            min_gap: gap,
            max_gap: gap + 1,
            phantom: PhantomData,
        }
    }

    /// generate all GappedKmerCtr's representing a fixed-length motif, given contant MIN_KMER_LEN
    /// eg (assuming MIN_KMER_LEN of 4), given motif_len 10, return motifs representing:
    /// (4, 0, 6), (4, 1, 5), (4, 2, 4), (5, 0, 5), (5, 1, 4), ...
    pub fn rect_tuples_for_motif_len(motif_len: usize) -> Result<Vec<(usize, usize, usize)>, ()> {
        if motif_len < 2 * MIN_KMER_LEN {
            return Err(());
        }

        let mut ctrs = vec![];
        debug!(
            "{:?}",
            MIN_KMER_LEN..1 + min(TOT_KMER_LEN, motif_len - MIN_KMER_LEN)
        );

        for gap_start in MIN_KMER_LEN..1 + min(motif_len, TOT_KMER_LEN - MIN_KMER_LEN) {
            debug!(
                "{:?}",
                max(gap_start, motif_len - (TOT_KMER_LEN - gap_start))
                    ..max(gap_start, motif_len - MIN_KMER_LEN)
            );

            for gap_end in max(gap_start, motif_len - (TOT_KMER_LEN - gap_start))
                ..max(gap_start, motif_len - MIN_KMER_LEN)
            {
                ctrs.push((gap_start, gap_end - gap_start, motif_len - gap_end));
            }
        }
        Ok(ctrs)
    }

    /// retrieve value at position, where kmers are specified as slices.  gap is in sequence space, not underlying array space.
    /// eg, foo.get(b"ATGC", b"TTCA", 0) -> 10
    pub fn get(&self, kmer1: &[u8], kmer2: &[u8], gap: usize) -> usize {
        self.ctr[[
            GappedKmerCtr::<M>::kmer_to_int(kmer1),
            GappedKmerCtr::<M>::kmer_to_int(kmer2),
            gap - self.min_gap,
        ]]
    }

    /// increment cell and return new value; indices are specified as they are for get
    pub fn incr(&mut self, kmer1: &[u8], kmer2: &[u8], gap: usize) -> usize {
        let idx = [
            GappedKmerCtr::<DNAMotif>::kmer_to_int(kmer1),
            GappedKmerCtr::<DNAMotif>::kmer_to_int(kmer2),
            gap - self.min_gap,
        ];
        self.ctr[idx] += 1;
        self.ctr[idx]
    }

    /// given a sequence, increment all kmers
    pub fn update_with_seq(&mut self, seq: &[u8]) {
        /*
        // annoying - work around borrow-checking
        let kmer_len = self.kmer_len;
        let min_gap = self.min_gap;
        let max_gap = self.max_gap;
         */
        for gap in self.min_gap..self.max_gap {
            for start_kmer1 in 0..1 + seq.len() - (self.kmer_len.0 + self.kmer_len.1 + gap) {
                let end_kmer1 = start_kmer1 + self.kmer_len.0;
                self.incr(
                    &seq[start_kmer1..end_kmer1],
                    &seq[end_kmer1 + gap..end_kmer1 + gap + self.kmer_len.1],
                    gap,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEQ: &'static [u8] = b"TGCACGGT";
    #[test]
    fn kmer_to_int() {
        assert_eq!(GappedKmerCtr::<DNAMotif>::kmer_to_int(b"A"), 0);
        assert_eq!(GappedKmerCtr::<DNAMotif>::kmer_to_int(b"T"), 1);
        assert_eq!(GappedKmerCtr::<DNAMotif>::kmer_to_int(b"G"), 2);
        assert_eq!(GappedKmerCtr::<DNAMotif>::kmer_to_int(b"C"), 3);
        assert_eq!(GappedKmerCtr::<DNAMotif>::kmer_to_int(b"AA"), 0);
        assert_eq!(GappedKmerCtr::<DNAMotif>::kmer_to_int(b"TA"), 4);

        assert_eq!(GappedKmerCtr::<DNAMotif>::int_to_kmer(1, 0), b"A");
        assert_eq!(GappedKmerCtr::<DNAMotif>::int_to_kmer(1, 1), b"T");
        assert_eq!(GappedKmerCtr::<DNAMotif>::int_to_kmer(1, 2), b"G");
        assert_eq!(GappedKmerCtr::<DNAMotif>::int_to_kmer(1, 3), b"C");
        assert_eq!(GappedKmerCtr::<DNAMotif>::int_to_kmer(2, 0), b"AA");
        assert_eq!(GappedKmerCtr::<DNAMotif>::int_to_kmer(2, 4), b"TA");
    }

    #[test]
    fn index() {
        // TGC -> 27
        // CGG -> 58
        // GCA -> 44
        // GGT -> 41

        let mut ctr = GappedKmerCtr::<DNAMotif>::new(3, 1, 2);
        ctr.update_with_seq(SEQ);

        assert_eq!(ctr.get(b"TGC", b"CGG", 1), 1);

        for i in 0..64 {
            for j in 0..64 {
                assert_eq!(
                    ctr.ctr[[i, j, 0]],
                    match (i, j) {
                        (27, 58) => 1,
                        (44, 41) => 1,
                        _ => 0,
                    }
                );
            }
        }
    }

    #[test]
    fn test_rect() {
        // TGCACGGT
        let mut ctr = GappedKmerCtr::<DNAMotif>::rect(3, 1, 4);
        ctr.update_with_seq(SEQ);
        assert_eq!(ctr.get(b"TGC", b"CGGT", 1), 1);
    }

    #[test]
    fn test_rects_for_motif_len() {
        assert!(
            GappedKmerCtr::<DNAMotif>::rect_tuples_for_motif_len(2 * MIN_KMER_LEN - 1).is_err()
        );

        assert_eq!(
            GappedKmerCtr::<DNAMotif>::rect_tuples_for_motif_len(2 * MIN_KMER_LEN)
                .unwrap()
                .len(),
            1
        );
        let v = GappedKmerCtr::<DNAMotif>::rect_tuples_for_motif_len(2 * MIN_KMER_LEN + 1).unwrap();
        assert!(v.iter().any(|t| *t == (MIN_KMER_LEN, 0, MIN_KMER_LEN)));
    }
}

/*

impl AsRef<[usize; 3]> for &(&[u8], &[u8], usize) {
    fn as_ref(&self) -> &T;
}

impl Index<KmerCtrIdx> for GappedKmerCtr {
    type Output = usize;

    fn index(&self, idx: KmerCtrIdx) -> &usize {
        self.ctr[[0,0,0]]
    }
}

impl From<
*/
