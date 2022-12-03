use std::ops::Range;
use bio::data_structures::interval_tree::*;
use std::collections::{BTreeMap, HashMap};
use num::Integer;
use num_traits::{Num, NumAssignOps, NumCast};

use super::{GenomicRange, BEDLike, split_by_len};

pub struct BedTree<D>(HashMap<String, IntervalTree<u64, D>>);

impl<D> Default for BedTree<D> {
    fn default() -> Self {
        Self(HashMap::new())
    }
}

impl<D, B: BEDLike> FromIterator<(B, D)> for BedTree<D> {
    fn from_iter<I: IntoIterator<Item = (B, D)>>(iter: I) -> Self {
        let mut hmap: HashMap<String, Vec<(Range<u64>, D)>> = HashMap::new();
        for (bed, data) in iter {
            let chr = bed.chrom();
            let interval = bed.start() .. bed.end();
            let vec = hmap.entry(chr.to_string()).or_insert(Vec::new());
            vec.push((interval, data));
        }
        //let hm = hmap.into_iter().map(|(chr, vec)| (chr, vec.into_iter().unique().collect())).collect();
        let hm = hmap.into_iter().map(|(chr, vec)| (chr, vec.into_iter().collect())).collect();
        BedTree(hm)
    }
}

/// An `IntervalTreeIterator` is returned by `Intervaltree::find` and iterates over the entries
/// overlapping the query
pub struct BedTreeIterator<'a, D> {
    chrom: String,
    interval_tree_iterator: Option<IntervalTreeIterator<'a, u64, D>>,
}

impl<'a, D: 'a> Iterator for BedTreeIterator<'a, D> {
    type Item = (GenomicRange, &'a D);

    fn next(&mut self) -> Option<Self::Item> {
        match self.interval_tree_iterator {
            None => return None,
            Some(ref mut iter) => match iter.next() {
                None => return None,
                Some(item) => {
                    let bed = GenomicRange::new(self.chrom.to_string(), item.interval().start, item.interval().end);
                    Some((bed, item.data()))
                }
            }
        }
    }
}


impl<D> BedTree<D> {
    pub fn find<B: BEDLike>(&self, bed: &B) -> BedTreeIterator<'_, D> {
        let chr = bed.chrom().to_string();
        let interval = bed.start() .. bed.end();
        match self.0.get(&chr) {
            None => BedTreeIterator { chrom: chr, interval_tree_iterator: None },
            Some(tree) => BedTreeIterator { chrom: chr, interval_tree_iterator: Some(tree.find(interval)) }
        }
    }

    pub fn is_overlapped<B: BEDLike>(&self, bed: &B) -> bool {
        self.find(bed).next().is_some()
    }
}

pub struct GenomeRegions<B> {
    pub regions: Vec<B>,
    pub indices: BedTree<usize>,
}

impl<B: BEDLike> GenomeRegions<B> {
    pub fn get_regions(&self) -> &Vec<B> { &self.regions }

    pub fn len(&self) -> usize { self.regions.len() }
}

impl<B: BEDLike> FromIterator<B> for GenomeRegions<B> {
    fn from_iter<I: IntoIterator<Item = B>>(iter: I) -> Self {
        let regions: Vec<B> = iter.into_iter().collect();
        let indices = regions.iter().enumerate().map(|(i, b)| (b.to_genomic_range(), i)).collect();
        GenomeRegions { regions, indices }
    }
}

#[derive(Clone)]
pub struct Coverage<'a, B, N> {
    pub consumed_tags: f64,
    genome_regions: &'a GenomeRegions<B>,
    coverage: Vec<N>,
}

impl <'a, N: Num + NumCast + NumAssignOps + Copy, B: BEDLike> Coverage<'a, B, N> {
    pub fn new(genome_regions: &'a GenomeRegions<B>) -> Self {
        Self {
            consumed_tags: 0.0,
            genome_regions,
            coverage: vec![N::zero(); genome_regions.len()],
        }
    }

    pub fn reset(&mut self) {
        self.consumed_tags = 0.0;
        self.coverage.fill(N::zero());
    }

    pub fn insert<D>(&mut self, tag: &D, count: N)
    where
        D: BEDLike,
    {
        self.consumed_tags += <f64 as NumCast>::from(count).unwrap();
        self.genome_regions.indices.find(tag).for_each(|(_, idx)| self.coverage[*idx] += count);
    }

    pub fn get_regions(&'a self) -> impl Iterator<Item = &B> + 'a
    {
        self.genome_regions.regions.iter()
    }

    pub fn get_coverage(&self) -> &Vec<N> { &self.coverage }
}


#[derive(Clone)]
pub struct SparseCoverage<'a, B, N> {
    pub consumed_tags: f64,
    genome_regions: &'a GenomeRegions<B>,
    coverage: BTreeMap<usize, N>,
}

impl <'a, N: Num + NumCast + NumAssignOps + Copy, B: BEDLike> SparseCoverage<'a, B, N> {
    pub fn new(genome_regions: &'a GenomeRegions<B>) -> Self {
        Self {
            consumed_tags: 0.0,
            genome_regions,
            coverage: BTreeMap::new(),
        }
    }

    pub fn reset(&mut self) {
        self.consumed_tags = 0.0;
        self.coverage = BTreeMap::new();
    }

    pub fn insert<D>(&mut self, tag: &D, count: N)
    where
        D: BEDLike,
    {
        self.consumed_tags += <f64 as NumCast>::from(count).unwrap();
        self.genome_regions.indices.find(tag).for_each(|(_, idx)| {
            let counter = self.coverage.entry(*idx).or_insert(N::zero());
            *counter += count;
        });
    }

    pub fn get_regions(&'a self) -> impl Iterator<Item = &B> + 'a
    {
        self.genome_regions.regions.iter()
    }

    pub fn get_coverage(&self) -> &BTreeMap<usize, N> { &self.coverage }

    pub fn get_coverage_as_vec(&self) -> Vec<N> {
        let mut coverage = vec![N::zero(); self.genome_regions.len()];
        self.coverage.iter().for_each(|(idx, v)| coverage[*idx] = *v);
        coverage
    }
}

#[derive(Clone)]
pub struct BinnedCoverage<'a, B, N> {
    pub len: usize,
    pub bin_size: u64,
    pub consumed_tags: f64,
    genome_regions: &'a GenomeRegions<B>,
    coverage: Vec<Vec<N>>,
}

impl <'a, N: Num + NumCast + NumAssignOps + Copy, B: BEDLike> BinnedCoverage<'a, B, N> {
    pub fn new(genome_regions: &'a GenomeRegions<B>, bin_size: u64) -> Self {
        let coverage: Vec<Vec<N>> = genome_regions.regions.iter()
            .map(|x| vec![N::zero(); x.len().div_ceil(&bin_size) as usize])
            .collect();
        let len = coverage.iter().map(|x| x.len()).sum();
        Self {genome_regions, len, bin_size, coverage, consumed_tags: 0.0}
    }

    pub fn reset(&mut self) {
        self.consumed_tags = 0.0;
        self.coverage.iter_mut().for_each(|x| x.fill(N::zero()));
    }

    pub fn insert<D>(&mut self, tag: &D, count: N)
    where
        D: BEDLike,
    {
        self.consumed_tags += <f64 as NumCast>::from(count).unwrap();
        self.genome_regions.indices.find(tag).for_each(|(region, out_idx)| {
            let i = tag.start().saturating_sub(region.start()).div_floor(&self.bin_size);
            let j = (tag.end() - 1 - region.start())
                .min(region.len() - 1).div_floor(&self.bin_size);
            (i..=j).for_each(|in_idx| self.coverage[*out_idx][in_idx as usize] += count);
        });
    }

    pub fn get_regions(&'a self) -> impl Iterator<Item = impl Iterator<Item = B>> + 'a
    where
        B: Clone,
    {
        self.genome_regions.regions.iter()
            .map(|x| super::split_by_len(x, self.bin_size))
    }

    pub fn get_coverage(&self) -> &Vec<Vec<N>> { &self.coverage }
}

#[derive(Clone)]
pub struct SparseBinnedCoverage<'a, B, N> {
    pub len: usize,
    pub bin_size: u64,
    pub consumed_tags: f64,
    genome_regions: &'a GenomeRegions<B>,
    accu_size: Vec<usize>,
    coverage: BTreeMap<usize, N>,
}

impl <'a, N: Num + NumCast + NumAssignOps + Copy, B: BEDLike> SparseBinnedCoverage<'a, B, N> {
    pub fn new(genome_regions: &'a GenomeRegions<B>, bin_size: u64) -> Self {
        let mut len = 0;
        let accu_size = genome_regions.regions.iter().map(|x| {
            let n = x.len().div_ceil(&bin_size) as usize;
            let output = len;
            len += n;
            output
        }).collect();
        Self {
            len, bin_size, consumed_tags: 0.0, genome_regions, accu_size,
            coverage: BTreeMap::new()
        }
    }

    pub fn reset(&mut self) {
        self.consumed_tags = 0.0;
        self.coverage = BTreeMap::new();
    }

    pub fn insert<D>(&mut self, tag: &D, count: N)
    where
        D: BEDLike,
    {
        self.consumed_tags += <f64 as NumCast>::from(count).unwrap();
        self.genome_regions.indices.find(tag).for_each(|(region, out_idx)| {
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

    pub fn get_regions(&'a self) -> impl Iterator<Item = impl Iterator<Item = B>> + 'a
    where
        B: Clone,
    {
        self.genome_regions.regions.iter()
            .map(|x| split_by_len(x, self.bin_size))
    }

    pub fn get_coverage(&self) -> &BTreeMap<usize, N> { &self.coverage }

    pub fn get_coverage_as_vec(&self) -> Vec<N> {
        let mut coverage = vec![N::zero(); self.len];
        self.coverage.iter().for_each(|(idx, v)| coverage[*idx] = *v);
        coverage
    }
}

#[cfg(test)]
mod bed_intersect_tests {
    use super::*;

    #[test]
    fn test_intersect() {
        let bed_set1: Vec<GenomicRange> = [
            "chr1:200-500",
            "chr1:200-500",
            "chr1:1000-2000",
        ].into_iter().map(|x| x.parse().unwrap()).collect();
        let bed_set2: Vec<GenomicRange> = [
            "chr1:100-210",
            "chr1:100-200",
        ].into_iter().map(|x| x.parse().unwrap()).collect();
        let tree: BedTree<()> = bed_set1.clone().into_iter().map(|x| (x, ())).collect();

        assert_eq!(
            tree.find(&bed_set2[0]).map(|x| x.0).collect::<Vec<_>>(),
            vec![bed_set1.as_slice()[0].clone(), bed_set1.as_slice()[0].clone()]
        );

        let result: Vec<GenomicRange> = bed_set2.into_iter().filter(|x| tree.is_overlapped(x)).collect();
        let expected = vec!["chr1:100-210".parse().unwrap()];
        assert_eq!(result, expected);

    }

    #[test]
    fn test_coverage() {
        let regions: GenomeRegions<GenomicRange> = vec![
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
    }

}