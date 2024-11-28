use itertools::Itertools;
use num::Num;
use num_traits::NumAssignOps;

use super::BedGraph;
use crate::bed::{GenomicRange, Score, Strand};

use std::{cmp::Ordering, iter::Sum, ops::Neg};

/// Common BED fields
pub trait BEDLike {
    /// Return the chromosome name of the record
    fn chrom(&self) -> &str;

    /// Change the chromosome name of the record
    fn set_chrom(&mut self, chrom: &str) -> &mut Self;

    /// Return the 0-based start position of the record
    fn start(&self) -> u64;

    /// Change the 0-based start position of the record
    fn set_start(&mut self, start: u64) -> &mut Self;

    /// Return the end position (non-inclusive) of the record
    fn end(&self) -> u64;

    /// Change the end position (non-inclusive) of the record
    fn set_end(&mut self, end: u64) -> &mut Self;

    /// Return the name of the record
    fn name(&self) -> Option<&str> {
        None
    }

    /// Return the score of the record
    fn score(&self) -> Option<Score> {
        None
    }

    /// Return the strand of the record
    fn strand(&self) -> Option<Strand> {
        None
    }

    /// Return the length of the record. Return 0 if the end position is smaller
    /// than the start position.
    fn len(&self) -> u64 {
        self.end().saturating_sub(self.start())
    }

    fn compare(&self, other: &Self) -> Ordering {
        self.chrom()
            .cmp(other.chrom())
            .then(self.start().cmp(&other.start()))
            .then(self.end().cmp(&other.end()))
    }

    /// Return the overlap
    fn overlap<B: BEDLike>(&self, other: &B) -> Option<GenomicRange> {
        if self.chrom() != other.chrom() {
            None
        } else {
            let start = self.start().max(other.start());
            let end = self.end().min(other.end());
            if start >= end {
                None
            } else {
                Some(GenomicRange::new(self.chrom(), start, end))
            }
        }
    }

    /// Return the size of overlap between two records
    fn n_overlap<B: BEDLike>(&self, other: &B) -> u64 {
        self.overlap(other).map_or(0, |x| x.len())
    }

    /// Convert the record to a `GenomicRange`
    fn to_genomic_range(&self) -> GenomicRange {
        GenomicRange::new(self.chrom(), self.start(), self.end())
    }

    /// Split into consecutive records with the specified length. The length of
    /// the last record may be shorter.
    fn split_by_len(&self, bin_size: u64) -> impl Iterator<Item = GenomicRange> {
        let start = self.start();
        let end = self.end();
        (start..end)
            .step_by(bin_size as usize)
            .map(move |x| GenomicRange::new(self.chrom(), x, (x + bin_size).min(end)))
    }

    /// Split into consecutive records with the specified length starting from the end.
    /// The result is in reverse order compared to `split_by_len`. The length of the last
    /// record may be shorter.
    fn rsplit_by_len(&self, bin_size: u64) -> impl Iterator<Item = GenomicRange> {
        let start = self.start();
        let end = self.end();
        (start + 1..=end)
            .rev()
            .step_by(bin_size as usize)
            .map(move |x| GenomicRange::new(self.chrom(), x.saturating_sub(bin_size).max(start), x))
    }
}

pub struct MergeBed<I, B, F> {
    sorted_bed_iter: I,
    merger: F,
    accum: Option<((String, u64, u64), Vec<B>)>,
}

impl<I, F, B, O> Iterator for MergeBed<I, B, F>
where
    I: Iterator<Item = B>,
    B: BEDLike,
    F: Fn(Vec<B>) -> O,
{
    type Item = O;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.sorted_bed_iter.next() {
                None => {
                    return if self.accum.is_none() {
                        None
                    } else {
                        let (_, accum) = std::mem::replace(&mut self.accum, None).unwrap();
                        Some((self.merger)(accum))
                    }
                }
                Some(record) => match self.accum.as_mut() {
                    None => {
                        self.accum = Some((
                            (record.chrom().to_string(), record.start(), record.end()),
                            vec![record],
                        ))
                    }
                    Some(((chr, s, e), accum)) => {
                        let chr_ = record.chrom();
                        let s_ = record.start();
                        let e_ = record.end();
                        if chr != chr_ || s_ > *e {
                            let (_, acc) = std::mem::replace(
                                &mut self.accum,
                                Some(((chr_.to_string(), s_, e_), vec![record])),
                            )
                            .unwrap();
                            return Some((self.merger)(acc));
                        } else if s_ < *s {
                            panic!("input is not sorted")
                        } else if e_ > *e {
                            *e = e_;
                            accum.push(record);
                        } else {
                            accum.push(record);
                        }
                    }
                },
            }
        }
    }
}

/// Merge sorted BED records. Overlapping records are processed according to the
/// function provided.
pub fn merge_sorted_bed_with<In, I, B, O, F>(sorted_bed_iter: In, merger: F) -> MergeBed<I, B, F>
where
    In: IntoIterator<IntoIter = I>,
    I: Iterator<Item = B>,
    B: BEDLike,
    F: Fn(Vec<B>) -> O,
{
    MergeBed {
        sorted_bed_iter: sorted_bed_iter.into_iter(),
        merger,
        accum: None,
    }
}

/// Merge sorted BED records. Overlapping records are concatenated into a single
/// record.
pub fn merge_sorted_bed<I, B>(sorted_iter: I) -> impl Iterator<Item = GenomicRange>
where
    I: IntoIterator<Item = B>,
    B: BEDLike,
{
    merge_sorted_bed_with(sorted_iter, |x| {
        GenomicRange::new(
            x[0].chrom(),
            x.iter().map(|x| x.start()).min().unwrap(),
            x.iter().map(|x| x.end()).max().unwrap(),
        )
    })
}

pub fn merge_sorted_bedgraph<I, V>(sorted_iter: I) -> impl Iterator<Item = BedGraph<V>>
where
    I: IntoIterator<Item = BedGraph<V>>,
    V: Num + NumAssignOps + Sum + Neg<Output = V> + PartialOrd + Copy,
{
    merge_sorted_bed_with(sorted_iter, |bdgs| {
        // group by start and end points
        let point_groups = bdgs
            .iter()
            .flat_map(|bedgraph| {
                [
                    (bedgraph.start(), bedgraph.value),
                    (bedgraph.end(), bedgraph.value.neg()),
                ]
            })
            .sorted_unstable_by_key(|x| x.0)
            .chunk_by(|x| (x.0, x.1 < V::zero()));
        let mut point_groups = point_groups.into_iter();

        let chrom = bdgs[0].chrom();
        let ((mut prev_pos, _), first_group) = point_groups.next().unwrap();
        let mut acc_val = first_group.into_iter().map(|x| x.1).sum();
        let mut prev_bedgraph = BedGraph::new(chrom, prev_pos, prev_pos, acc_val);

        let mut result = point_groups
            .flat_map(|((pos, _), group)| {
                let value_sum = group.into_iter().map(|x| x.1).sum();
                let mut bedgraph = None;

                if prev_pos != pos {
                    if acc_val == prev_bedgraph.value {
                        prev_bedgraph.set_end(pos);
                    } else {
                        bedgraph = Some(prev_bedgraph.clone());
                        prev_bedgraph = BedGraph::new(chrom, prev_pos, pos, acc_val);
                    }
                }

                acc_val += value_sum;
                prev_pos = pos;
                bedgraph
            })
            .collect::<Vec<_>>();
        result.push(prev_bedgraph);
        result
    })
    .flatten()
}

#[cfg(test)]
mod bed_tests {
    use super::*;
    use crate::extsort::ExternalSorterBuilder;
    use itertools::Itertools;
    use rand::Rng;

    fn rand_bed(n: usize) -> Vec<GenomicRange> {
        let mut rng = rand::thread_rng();
        let mut result = Vec::with_capacity(n);
        for _ in 0..n {
            let bed = GenomicRange::new(
                format!("chr{}", rng.gen_range(1..20)),
                rng.gen_range(0..10000),
                rng.gen_range(1000000..2000000),
            );
            result.push(bed);
        }
        result
    }

    #[test]
    fn test_sort() {
        let data1 = vec![
            GenomicRange::new("chr1", 1020, 1230),
            GenomicRange::new("chr1", 0, 500),
            GenomicRange::new("chr2", 500, 1000),
            GenomicRange::new("chr1", 500, 1000),
            GenomicRange::new("chr1", 1000, 1230),
            GenomicRange::new("chr1", 1000, 1230),
        ];
        let data2 = rand_bed(100000);

        [data1, data2].into_iter().for_each(|mut data| {
            let sorted1 = ExternalSorterBuilder::new()
                .num_threads(2)
                .with_chunk_size(10000)
                .build()
                .unwrap()
                .sort(data.clone().into_iter())
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            let sorted2 = ExternalSorterBuilder::new()
                .with_chunk_size(10000)
                .with_compression(4)
                .build()
                .unwrap()
                .sort(data.clone().into_iter())
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();
            data.sort();
            assert_eq!(data, sorted1);
            assert_eq!(data, sorted2);
        })
    }

    #[test]
    fn test_split() {
        assert_eq!(
            GenomicRange::new("chr1", 0, 1230)
                .split_by_len(500)
                .collect::<Vec<_>>(),
            vec![
                GenomicRange::new("chr1", 0, 500),
                GenomicRange::new("chr1", 500, 1000),
                GenomicRange::new("chr1", 1000, 1230),
            ],
        );
        assert_eq!(
            GenomicRange::new("chr1", 0, 1230)
                .rsplit_by_len(500)
                .collect::<Vec<_>>(),
            vec![
                GenomicRange::new("chr1", 730, 1230),
                GenomicRange::new("chr1", 230, 730),
                GenomicRange::new("chr1", 0, 230),
            ],
        );
        assert_eq!(
            GenomicRange::new("chr1", 500, 2500)
                .split_by_len(500)
                .collect::<Vec<_>>(),
            GenomicRange::new("chr1", 500, 2500)
                .rsplit_by_len(500)
                .sorted()
                .collect::<Vec<_>>(),
        );
        assert_eq!(
            GenomicRange::new("chr1", 500, 2500)
                .split_by_len(1)
                .collect::<Vec<_>>(),
            GenomicRange::new("chr1", 500, 2500)
                .rsplit_by_len(1)
                .sorted()
                .collect::<Vec<_>>(),
        );
    }

    #[test]
    fn test_overlap() {
        assert_eq!(
            GenomicRange::new("chr1", 100, 200).overlap(&GenomicRange::new("chr1", 150, 300)),
            Some(GenomicRange::new("chr1", 150, 200)),
        );

        assert_eq!(
            GenomicRange::new("chr1", 100, 200).overlap(&GenomicRange::new("chr1", 110, 190)),
            Some(GenomicRange::new("chr1", 110, 190)),
        );

        assert_eq!(
            GenomicRange::new("chr1", 100, 200).overlap(&GenomicRange::new("chr1", 90, 99)),
            None,
        );

        assert_eq!(
            GenomicRange::new("chr1", 111, 180).overlap(&GenomicRange::new("chr1", 110, 190)),
            Some(GenomicRange::new("chr1", 111, 180)),
        );

        assert_eq!(
            GenomicRange::new("chr1", 111, 200).overlap(&GenomicRange::new("chr1", 110, 190)),
            Some(GenomicRange::new("chr1", 111, 190)),
        );
    }

    #[test]
    fn test_merge() {
        let input = [
            (0, 100),
            (10, 20),
            (50, 150),
            (120, 160),
            (155, 200),
            (155, 220),
            (500, 1000),
            (2000, 2100),
            (2100, 2200),
        ]
        .into_iter()
        .map(|(a, b)| std::io::Result::Ok(GenomicRange::new("chr1", a, b)));
        let expect: Vec<GenomicRange> = [(0, 220), (500, 1000), (2000, 2200)]
            .into_iter()
            .map(|(a, b)| GenomicRange::new("chr1", a, b))
            .collect();
        assert_eq!(
            merge_sorted_bed(input.clone().map(|x| x.unwrap())).collect::<Vec<_>>(),
            expect
        );
    }

    #[test]
    fn test_merge_bedgraph() {
        let reads = vec![
            BedGraph::new("chr1", 0, 100, 1),
            BedGraph::new("chr1", 0, 100, 1),
            BedGraph::new("chr1", 2, 100, 2),
            BedGraph::new("chr1", 30, 60, 3),
            BedGraph::new("chr1", 40, 50, 4),
            BedGraph::new("chr1", 60, 65, 5),
        ];

        let actual: Vec<_> = merge_sorted_bedgraph(reads).collect();
        let expected = vec![
            BedGraph::new("chr1", 0, 2, 2),
            BedGraph::new("chr1", 2, 30, 4),
            BedGraph::new("chr1", 30, 40, 7),
            BedGraph::new("chr1", 40, 50, 11),
            BedGraph::new("chr1", 50, 60, 7),
            BedGraph::new("chr1", 60, 65, 9),
            BedGraph::new("chr1", 65, 100, 4),
        ];
        assert_eq!(expected, actual);
    }
}
