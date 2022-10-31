use crate::bed::{GenomicRange, Score, Strand};

use itertools::Itertools;
use std::{cmp::Ordering, io::Error};
use extsort::{ExternalSorter, sorter::Sortable};
use std::path::PathBuf;

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
    fn name(&self) -> Option<&str> { None }

    /// Return the score of the record
    fn score(&self) -> Option<Score> { None }

    /// Return the strand of the record
    fn strand(&self) -> Option<Strand> { None }

    /// Return the length of the record
    fn len(&self) -> u64 { self.end() - self.start() }

    fn compare(&self, other: &Self) -> Ordering {
        self.chrom().cmp(other.chrom())
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
}

/// Split into consecutive records with the specified length. The length of
/// the last record may be shorter.
pub fn split_by_len<B>(bed: &B, bin_size: u64) -> impl Iterator<Item = B>
where
    B: BEDLike + Clone,
{
    let start = bed.start();
    let end = bed.end();
    let mut bed_ = (*bed).clone();
    (start .. end).step_by(bin_size as usize).map(move |x| {
        bed_.set_start(x).set_end((x + bin_size).min(end));
        bed_.clone()
    })
}

pub fn sort_bed<I, B>(bed_iter: I) -> impl Iterator<Item = B>
where
    I: Iterator<Item = B>,
    B: BEDLike + Sortable,
{
    ExternalSorter::new()
        .with_segment_size(50000000)
        .with_sort_dir(PathBuf::from("./"))
        .with_parallel_sort()
        .sort_by(bed_iter, BEDLike::compare).unwrap()
}

pub fn sort_bed_by<I, B, F>(bed_iter: I, cmp: F) -> impl Iterator<Item = B>
where
    I: Iterator<Item = B>,
    B: BEDLike + Sortable,
    F: Fn(&B, &B) -> Ordering + Send + Sync, 
{
    ExternalSorter::new()
        .with_segment_size(50000000)
        .with_sort_dir(PathBuf::from("./"))
        .with_parallel_sort()
        .sort_by(bed_iter, cmp).unwrap()
}

pub fn sort_bed_by_key<I, B, F, K>(bed_iter: I, f: F) -> impl Iterator<Item = B>
where
    I: Iterator<Item = B>,
    B: BEDLike + Sortable,
    F: Fn(&B) -> K + Send + Sync,
    K: Ord,
{
    ExternalSorter::new()
        .with_segment_size(50000000)
        .with_sort_dir(PathBuf::from("./"))
        .with_parallel_sort()
        .sort_by_key(bed_iter, f).unwrap()
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
                None => return if self.accum.is_none() {
                    None
                } else {
                    let (_, accum) = std::mem::replace(&mut self.accum, None).unwrap();
                    Some((self.merger)(accum))
                },
                Some(record) => match self.accum.as_mut() {
                    None => self.accum = Some((
                        (record.chrom().to_string(), record.start(), record.end()),
                        vec![record],
                    )),
                    Some(((chr, s, e), accum)) => {
                        let chr_ = record.chrom();
                        let s_ = record.start();
                        let e_ = record.end();
                        if chr != chr_ || s_ > *e {
                            let (_, acc) = std::mem::replace(
                                &mut self.accum,
                                Some(((chr_.to_string(), s_, e_), vec![record])),
                            ).unwrap();
                            return Some((self.merger)(acc))
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

pub fn merge_sorted_bed_with<I, B, O, F>(sorted_bed_iter: I, merger: F) -> MergeBed<I, B, F>
where
    I: Iterator<Item = B>,
    B: BEDLike,
    F: Fn(Vec<B>) -> O,
{
    MergeBed { sorted_bed_iter, merger, accum: None }
}

pub fn merge_bed_with<I, B, O, F>(bed_iter: I, merger: F) -> MergeBed<impl Iterator<Item = B>, B, F>
where
    I: Iterator<Item = B>,
    B: BEDLike + Sortable,
    F: Fn(Vec<B>) -> O,
{
    merge_sorted_bed_with(sort_bed(bed_iter), merger)
}

pub fn merge_sorted_bed<I, B>(sorted_bed_iter: I) -> MergeBed<I, B, impl Fn(Vec<B>) -> GenomicRange>
where
    I: Iterator<Item = B>,
    B: BEDLike,
{
    let merger = |x: Vec<B>| -> GenomicRange { GenomicRange::new(
        x[0].chrom(),
        x.iter().map(|x| x.start()).min().unwrap(),
        x.iter().map(|x| x.end()).max().unwrap(),
    ) };
    MergeBed { sorted_bed_iter, merger, accum: None }
}

pub fn merge_bed<I, B>(
    bed_iter: I
) -> MergeBed<impl Iterator<Item = B>, B, impl Fn(Vec<B>) -> GenomicRange>
where
    I: Iterator<Item = B>,
    B: BEDLike + Sortable,
{
    merge_sorted_bed(sort_bed(bed_iter))
}


#[cfg(test)]
mod bed_tests {
    use super::*;

    #[test]
    fn test_sort() {
        let mut data = vec![
            GenomicRange::new("chr1", 1020, 1230),
            GenomicRange::new("chr1", 0, 500),
            GenomicRange::new("chr2", 500, 1000),
            GenomicRange::new("chr1", 500, 1000),
            GenomicRange::new("chr1", 1000, 1230),
            GenomicRange::new("chr1", 1000, 1230),
        ];
        let sorted: Vec<_> = sort_bed(data.clone().into_iter()).collect();
        data.sort();
        assert_eq!(data, sorted);
    }

    #[test]
    fn test_split() {
        let beds: Vec<GenomicRange> = split_by_len(&GenomicRange::new("chr1", 0, 1230), 500).collect();
        let expected = vec![
            GenomicRange::new("chr1", 0, 500),
            GenomicRange::new("chr1", 500, 1000),
            GenomicRange::new("chr1", 1000, 1230),
        ];
        assert_eq!(beds, expected);
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
        ].into_iter().map(|(a,b)| GenomicRange::new("chr1", a, b));
        let expect: Vec<GenomicRange> = [
            (0, 220),
            (500, 1000),
            (2000, 2200),
        ].into_iter().map(|(a,b)| GenomicRange::new("chr1", a, b)).collect();
        assert_eq!(merge_sorted_bed(input.clone()).collect::<Vec<_>>(), expect);
        assert_eq!(merge_bed(input.rev()).collect::<Vec<_>>(), expect);
    }
}