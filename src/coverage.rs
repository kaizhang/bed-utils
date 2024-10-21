use std::collections::BTreeMap;
use num::Integer;
use num_traits::{Num, NumAssignOps, NumCast};

use crate::bed::{GenomicRange, BEDLike, split_by_len, map::GIntervalIndexSet};

#[derive(Debug, Clone)]
pub struct Coverage<'a, N> {
    total_count: f64,
    intervals: &'a GIntervalIndexSet,
    coverage: Vec<N>,
}

impl <'a, N: Num + NumCast + NumAssignOps + Copy> Coverage<'a, N> {
    pub fn new(intervals: &'a GIntervalIndexSet) -> Self {
        Self {
            total_count: 0.0,
            intervals,
            coverage: vec![N::zero(); intervals.len()],
        }
    }

    pub fn len(&self) -> usize { self.coverage.len() }

    pub fn total_count(&self) -> f64 { self.total_count }

    pub fn reset(&mut self) {
        self.total_count = 0.0;
        self.coverage.fill(N::zero());
    }

    pub fn get_index<D>(&self, tag: &D) -> impl Iterator<Item = usize> + '_
    where
        D: BEDLike,
    {
        self.intervals.find_index_of(tag)
    }

    pub fn insert<D>(&mut self, tag: &D, multiplicity: N)
    where
        D: BEDLike,
    {
        self.total_count += <f64 as NumCast>::from(multiplicity).unwrap();
        self.intervals.find_index_of(tag).for_each(|idx| self.coverage[idx] += multiplicity);
    }

    pub fn insert_at_index<D>(&mut self, index: usize, count: N)
    {
        self.total_count += <f64 as NumCast>::from(count).unwrap();
        self.coverage[index] += count;
    }

    pub fn regions(&'a self) -> impl Iterator<Item = &GenomicRange> + 'a
    {
        self.intervals.iter()
    }

    pub fn get_coverage(&self) -> &Vec<N> { &self.coverage }
}

#[derive(Clone)]
pub struct SparseCoverage<'a, N> {
    total_count: f64,
    intervals: &'a GIntervalIndexSet,
    coverage: BTreeMap<usize, N>,
}

impl <'a, N: Num + NumCast + NumAssignOps + Copy> SparseCoverage<'a, N> {
    pub fn new(intervals: &'a GIntervalIndexSet) -> Self {
        Self {
            total_count: 0.0,
            intervals,
            coverage: BTreeMap::new(),
        }
    }

    pub fn len(&self) -> usize { self.intervals.len() }

    pub fn total_count(&self) -> f64 { self.total_count }

    pub fn reset(&mut self) {
        self.total_count = 0.0;
        self.coverage = BTreeMap::new();
    }

    pub fn get_index<D>(&self, tag: &D) -> impl Iterator<Item = usize> + '_
    where
        D: BEDLike,
    {
        self.intervals.find_index_of(tag)
    }

    pub fn insert<D>(&mut self, tag: &D, count: N)
    where
        D: BEDLike,
    {
        self.total_count += <f64 as NumCast>::from(count).unwrap();
        self.intervals.find_index_of(tag).for_each(|idx| {
            self.coverage.entry(idx).and_modify(|x| *x += count).or_insert(count);
        });
    }

    pub fn insert_at_index<D>(&mut self, index: usize, count: N)
    {
        self.total_count += <f64 as NumCast>::from(count).unwrap();
        self.coverage.entry(index).and_modify(|x| *x += count).or_insert(count);
    }

    pub fn regions(&'a self) -> impl Iterator<Item = &GenomicRange> + 'a
    {
        self.intervals.iter()
    }

    pub fn get_coverage(&self) -> &BTreeMap<usize, N> { &self.coverage }

    pub fn get_coverage_as_vec(&self) -> Vec<N> {
        let mut coverage = vec![N::zero(); self.intervals.len()];
        self.coverage.iter().for_each(|(idx, v)| coverage[*idx] = *v);
        coverage
    }
}

#[derive(Debug, Clone)]
pub struct BinnedCoverage<'a, N> {
    len: usize,
    bin_size: u64,
    consumed_tags: f64,
    intervals: &'a GIntervalIndexSet,
    coverage: Vec<Vec<N>>,
}

impl <'a, N: Num + NumCast + NumAssignOps + Copy> BinnedCoverage<'a, N> {
    pub fn new(intervals: &'a GIntervalIndexSet, bin_size: u64) -> Self {
        let coverage: Vec<Vec<N>> = intervals.iter()
            .map(|x| vec![N::zero(); x.len().div_ceil(bin_size) as usize])
            .collect();
        let len = coverage.iter().map(|x| x.len()).sum();
        Self {intervals, len, bin_size, coverage, consumed_tags: 0.0}
    }

    pub fn len(&self) -> usize { self.len }

    pub fn total_count(&self) -> f64 { self.consumed_tags }

    pub fn reset(&mut self) {
        self.consumed_tags = 0.0;
        self.coverage.iter_mut().for_each(|x| x.fill(N::zero()));
    }

    pub fn insert<D>(&mut self, tag: &D, count: N)
    where
        D: BEDLike,
    {
        self.consumed_tags += <f64 as NumCast>::from(count).unwrap();
        self.intervals.find_full(tag).for_each(|(region, out_idx)| {
            let i = tag.start().saturating_sub(region.start()).div_floor(&self.bin_size);
            let j = (tag.end() - 1 - region.start())
                .min(region.len() - 1).div_floor(&self.bin_size);
            (i..=j).for_each(|in_idx| self.coverage[*out_idx][in_idx as usize] += count);
        });
    }

    pub fn regions(&'a self) -> impl Iterator<Item = impl Iterator<Item = GenomicRange>> + 'a
    {
        self.intervals.iter()
            .map(|x| split_by_len(x, self.bin_size))
    }

    pub fn get_coverage(&self) -> &Vec<Vec<N>> { &self.coverage }
}

#[derive(Debug, Clone)]
pub struct SparseBinnedCoverage<'a, N> {
    pub len: usize,
    pub bin_size: u64,
    pub consumed_tags: f64,
    intervals: &'a GIntervalIndexSet,
    accu_size: Vec<usize>,
    coverage: BTreeMap<usize, N>,
}

impl <'a, N: Num + NumCast + NumAssignOps + Copy> SparseBinnedCoverage<'a, N> {
    pub fn new(intervals: &'a GIntervalIndexSet, bin_size: u64) -> Self {
        let mut len = 0;
        let accu_size = intervals.iter().map(|x| {
            let n = x.len().div_ceil(bin_size) as usize;
            let output = len;
            len += n;
            output
        }).collect();
        Self {
            len, bin_size, consumed_tags: 0.0, intervals, accu_size,
            coverage: BTreeMap::new()
        }
    }

    pub fn len(&self) -> usize { self.len }

    pub fn total_count(&self) -> f64 { self.consumed_tags }

    pub fn reset(&mut self) {
        self.consumed_tags = 0.0;
        self.coverage = BTreeMap::new();
    }

    pub fn insert<D>(&mut self, tag: &D, count: N)
    where
        D: BEDLike,
    {
        self.consumed_tags += <f64 as NumCast>::from(count).unwrap();
        self.intervals.find_full(tag).for_each(|(region, out_idx)| {
            let i = tag.start().saturating_sub(region.start()).div_floor(&self.bin_size);
            let j = (tag.end() - 1 - region.start())
                .min(region.len() - 1).div_floor(&self.bin_size);
            let n = self.accu_size[*out_idx];
            (i..=j).for_each(|in_idx| {
                let counter = self.coverage.entry(n + in_idx as usize).or_insert(N::zero());
                *counter += count;
            });
        });
    }

    pub fn get_regions(&'a self) -> impl Iterator<Item = impl Iterator<Item = GenomicRange>> + 'a
    {
        self.intervals.iter()
            .map(|x| split_by_len(x, self.bin_size))
    }

    pub fn get_coverage(&self) -> &BTreeMap<usize, N> { &self.coverage }

    pub fn get_coverage_as_vec(&self) -> Vec<N> {
        let mut coverage = vec![N::zero(); self.len];
        self.coverage.iter().for_each(|(idx, v)| coverage[*idx] = *v);
        coverage
    }

    pub fn get_region_coverage(&'a self) -> impl Iterator<Item = (GenomicRange, N)> + 'a {
        self.get_coverage().iter().map(|(i, v)| {
            let region = match self.accu_size.binary_search(&i) {
                Ok(j) => {
                    let site = &self.intervals[j];
                    let chr = site.chrom();
                    let start = site.start();
                    let end = (start + self.bin_size).min(site.end());
                    GenomicRange::new(chr, start, end)
                },
                Err(j) => {
                    let site = &self.intervals[j-1];
                    let chr = site.chrom();
                    let prev = self.accu_size[j-1];
                    let start = site.start() + ((i - prev) as u64) * self.bin_size;
                    let end = (start + self.bin_size).min(site.end());
                    GenomicRange::new(chr, start, end)
                }
            };
            (region, *v)
        })
    }
}

#[cfg(test)]
mod bed_intersect_tests {
    use super::*;

    #[test]
    fn test_coverage() {
        let regions: GIntervalIndexSet = vec![
            GenomicRange::new("chr1".to_string(), 200, 500),
            GenomicRange::new("chr1".to_string(), 1000, 2000),
            GenomicRange::new("chr1".to_string(), 10000, 11000),
            GenomicRange::new("chr10".to_string(), 10, 20),
            GenomicRange::new("chr1".to_string(), 200, 500),
        ].into_iter().collect();
        let tags: Vec<GenomicRange> = vec![
            GenomicRange::new("chr1".to_string(), 100, 210),
            GenomicRange::new("chr1".to_string(), 100, 500),
            GenomicRange::new("chr1".to_string(), 100, 5000),
            GenomicRange::new("chr1".to_string(), 100, 200),
            GenomicRange::new("chr1".to_string(), 1000, 1001),
        ];

        let mut cov1 = Coverage::new(&regions);
        tags.iter().for_each(|x| cov1.insert(x, 1));
        let result1: Vec<u64> = cov1.get_coverage().to_vec();
        let mut cov2 = SparseCoverage::new(&regions);
        tags.iter().for_each(|x| cov2.insert(x, 1));
        let result2: Vec<u64> = cov2.get_coverage_as_vec();

        assert_eq!(result1, result2);
        assert_eq!(result1, vec![3, 2, 0, 0, 3]);

        let mut cov3 = BinnedCoverage::new(&regions, 100);
        tags.iter().for_each(|x| cov3.insert(x, 1));
        let result3: Vec<u64> = cov3.get_coverage().iter().flatten().map(|x| *x).collect();
        let mut cov4 = SparseBinnedCoverage::new(&regions, 100);
        tags.iter().for_each(|x| cov4.insert(x, 1));
        let result4: Vec<u64> = cov4.get_coverage_as_vec();

        assert_eq!(result3, result4);
        assert_eq!(
            result3,
            vec![
                3, 2, 2,
                2, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0,
                3, 2, 2,
            ]
        );
    }

    #[test]
    fn test_sparse_coverage() {
        let regions = [
            GenomicRange::new("chr1", 0, 2000),
            GenomicRange::new("chr2", 100, 2100),
            GenomicRange::new("chr3", 3000, 3500),
            ].into_iter().collect();
        let mut cov = SparseBinnedCoverage::new(&regions, 400);

        cov.insert(&GenomicRange::new("chr1", 3, 5), 1);
        cov.insert(&GenomicRange::new("chr2", 0, 500), 1);
        cov.insert(&GenomicRange::new("chr2", 0, 501), 1);
        cov.insert(&GenomicRange::new("chr3", 3400, 3401), 1);
        cov.insert(&GenomicRange::new("chr4", 0, 501), 1);
        cov.insert(&GenomicRange::new("chr3", 0, 501), 1);

        assert_eq!(
            vec![(0, 1), (5, 2), (6, 1), (11, 1)],
            cov.get_coverage().iter().map(|(i, x)| (*i, *x)).collect::<Vec<(usize, u64)>>()
        );

        let expected: Vec<_> = cov.get_regions().flatten().zip(
            cov.get_coverage_as_vec().iter()
        ).flat_map(|(region, x)| if *x == 0 {
            None
        } else {
            Some((region, *x))
        }).collect();
        assert_eq!(
            expected,
            cov.get_region_coverage().collect::<Vec<_>>(),
        );
    }
}