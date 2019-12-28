#[macro_use]
extern crate log;
extern crate bio;
extern crate gundam;

extern crate env_logger;
use bio::alignment::distance::hamming;
use bio::alignment::pairwise::*;
use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation::*;
use bio::utils::TextSlice;
use env_logger::Builder as LogBuilder;
use gundam::*;
use std::cmp::min;
use std::collections::HashSet;
use std::env;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::process::exit;
use std::vec;

#[derive(Debug)]
pub struct ArgError;
impl fmt::Display for ArgError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Incorrect number of arguments")
    }
}
impl Error for ArgError {}

type Cluster = Vec<Vec<u8>>;

fn ungapped_aln<F: MatchFunc>(
    x: TextSlice<'_>,
    y: TextSlice<'_>,
    aligner: &mut Aligner<F>,
    match_cutoff: f32,
    coverage_cutoff: f32,
) -> Option<(usize, usize, usize)> {
    let aln = aligner.local(x, y);
    let (offset_x, offset_y, len) = if aln.xstart > aln.ystart {
        let len = min(x.len() - (aln.xstart - aln.ystart), y.len());

        (aln.xstart - aln.ystart, 0, len)
    } else {
        let len = min(y.len() - (aln.ystart - aln.xstart), x.len());
        (0, aln.ystart - aln.xstart, len)
    };
    let fst = &x[offset_x..offset_x + len];
    let sec = &y[offset_y..offset_y + len];

    let match_ = (len as f32 - hamming(fst, sec) as f32) / len as f32;
    let coverage = (len as f32) / min(x.len(), y.len()) as f32;

    if match_ >= match_cutoff && coverage >= coverage_cutoff {
        Some((offset_x, offset_y, len))
    } else {
        None
    }
}

fn cluster<F: MatchFunc>(motifs: &Vec<Vec<u8>>, aligner: &mut Aligner<F>) -> Vec<Cluster> {
    let match_cutoff = 0.8;
    let coverage_cutoff = 0.8;
    let mut matches = vec![];
    for i in 0..motifs.len() - 1 {
        let v: Vec<_> = (i + 1..motifs.len())
            .filter_map(|j| {
                ungapped_aln(
                    &motifs[i],
                    &motifs[j],
                    aligner,
                    match_cutoff,
                    coverage_cutoff,
                )
                .and_then(|t| Some((j, t)))
            })
            .collect();
        matches.push(v);
    }

    if let Some((idx, ref biggest)) = matches.iter().enumerate().max_by_key(|(_, ref v)| v.len()) {
        if biggest.len() == 0 {
            // no more clusters left
            motifs
                .iter()
                .map(|ref m| vec![m.iter().cloned().collect()])
                .collect()
        } else {
            let claimed: HashSet<_> = biggest.iter().map(|(j, _)| j).collect();
            let mut this_clust: Cluster = vec![motifs[idx].iter().cloned().collect()];
            let mut rest = vec![];
            for (j, t) in biggest.into_iter().enumerate() {
                if claimed.contains(&j) {
                    this_clust.push(motifs[j].iter().cloned().collect());
                } else {
                    rest.push(motifs[j].iter().cloned().collect());
                }
            }
            let mut all = vec![this_clust];
            let other_clusters: Vec<Cluster> = cluster(&rest, aligner);
            all.extend(other_clusters.into_iter());
            all
        }
    } else {
        motifs
            .iter()
            .map(|ref m| vec![m.iter().cloned().collect()])
            .collect()
    }
}

fn main() -> Result<(), Box<Error>> {
    let _ = LogBuilder::new()
        .parse(&env::var("RUST_LOG").unwrap_or_default())
        .init();

    let args = env::args().collect::<Vec<String>>();
    if args.len() != 2 {
        println!(
            "invalid # of args: {}.  usage: gundam <motifs.txt>",
            args.len()
        );
        Err(ArgError)?
    }

    let file = File::open(&args[1]).expect("can't open index file");
    let idx_file = BufReader::new(&file);
    let motifs_: Result<Vec<Vec<u8>>, _> = idx_file
        .lines()
        .map(|line| line.and_then(|l| Ok(l.into_bytes())))
        .collect();
    let mut motifs = motifs_?;
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    let mut aligner = Aligner::new(-999, -99, &score);

    for clust in cluster(&motifs, &mut aligner) {
        for s in clust {
            println!("{}", Seq(s.as_ref()));
        }
        println!("###############\n");
    }

    Ok(())
}
