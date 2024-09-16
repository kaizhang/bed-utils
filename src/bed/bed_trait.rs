use crate::bed::{GenomicRange, Score, Strand};

use std::cmp::Ordering;

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

    /// Return the length of the record. Return 0 if the end position is smaller
    /// than the start position.
    fn len(&self) -> u64 { self.end().saturating_sub(self.start()) }

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

#[cfg(test)]
mod bed_tests {
    use crate::extsort::ExternalSorterBuilder;
    use rand::Rng;
    use super::*;

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
                .sort(data.clone().into_iter()).unwrap()
                .collect::<Result<Vec<_>, _>>().unwrap();
            let sorted2 = ExternalSorterBuilder::new()
                .with_chunk_size(10000)
                .with_compression(4)
                .build()
                .unwrap()
                .sort(data.clone().into_iter()).unwrap()
                .collect::<Result<Vec<_>, _>>().unwrap();
            data.sort();
            assert_eq!(data, sorted1);
            assert_eq!(data, sorted2);
        })
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
        ].into_iter().map(|(a,b)| std::io::Result::Ok(GenomicRange::new("chr1", a, b)));
        let expect: Vec<GenomicRange> = [
            (0, 220),
            (500, 1000),
            (2000, 2200),
        ].into_iter().map(|(a,b)| GenomicRange::new("chr1", a, b)).collect();
        assert_eq!(merge_sorted_bed(input.clone().map(|x| x.unwrap())).collect::<Vec<_>>(), expect);
    }
}