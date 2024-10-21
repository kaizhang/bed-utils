use std::ops::Index;
use std::collections::{HashMap, VecDeque};

use super::{GenomicRange, BEDLike};
use crate::intervaltree::{Interval, IterFind, IterLapper, Lapper};

/// A map from genomic internvals to values, internally represented as a interval tree
/// for efficient queries.
/// Note that GIntervalMap permits duplicated records.
#[derive(Debug, Clone)]
pub struct GIntervalMap<D>(HashMap<String, Lapper<u64, D>>);

impl<D> GIntervalMap<D> {
    pub fn new() -> Self {
        Self(HashMap::new())
    }

    /// Return the number of records.
    pub fn len(&self) -> usize {
        self.0.values().map(|x| x.len()).sum()
    }

    pub fn iter(&self) -> Iter<'_, D> {
        Iter(self.0.iter().map(|(k, v)| (k, v.iter())).collect::<VecDeque<_>>())
    }

    /// This is very inefficient and should be avoided if possible.
    pub fn insert<B: BEDLike>(&mut self, bed: &B, value: D) {
        let chr = bed.chrom();
        let interval = Interval { start: bed.start(), stop: bed.end(), val: value };
        let tree = self.0.entry(chr.to_string()).or_insert(Lapper::new(Vec::new()));
        tree.insert(interval);
    }

    /// Return regions that overlaop with the query.
    pub fn find<B: BEDLike>(&self, bed: &B) -> GIntervalQueryIter<'_, D> {
        let chr = bed.chrom().to_string();
        match self.0.get(&chr) {
            None => GIntervalQueryIter { chrom: chr, query_iter: None },
            Some(tree) => GIntervalQueryIter { chrom: chr, query_iter: Some(tree.find(bed.start(), bed.end())) }
        }
    }

    /// Determine if the query overlaps with any record.
    pub fn is_overlapped<B: BEDLike>(&self, bed: &B) -> bool {
        self.find(bed).next().is_some()
    }
}

impl<B: BEDLike, D> FromIterator<(B, D)> for GIntervalMap<D> {
    fn from_iter<I: IntoIterator<Item = (B, D)>>(iter: I) -> Self {
        let mut hmap: HashMap<String, Vec<_>> = HashMap::new();
        for (bed, data) in iter {
            let chr = bed.chrom();
            let interval = Interval { start: bed.start(), stop: bed.end(), val: data};
            let vec = hmap.entry(chr.to_string()).or_insert(Vec::new());
            vec.push(interval);
        }
        let hm = hmap.into_iter().map(|(chr, vec)| (chr, Lapper::new(vec))).collect();
        Self(hm)
    }
}

/// An `GIntervalQueryIter` is returned by `Intervaltree::find` and iterates over the entries
/// overlapping the query
#[derive(Debug)]
pub struct GIntervalQueryIter<'a, D> {
    chrom: String,
    query_iter: Option<IterFind<'a, u64, D>>,
}

impl<'a, D: 'a> Iterator for GIntervalQueryIter<'a, D> {
    type Item = (GenomicRange, &'a D);

    fn next(&mut self) -> Option<Self::Item> {
        match self.query_iter {
            None => return None,
            Some(ref mut iter) => match iter.next() {
                None => return None,
                Some(interval) => {
                    let bed = GenomicRange::new(&self.chrom, interval.start, interval.stop);
                    Some((bed, &interval.val))
                }
            }
        }
    }
}

pub struct Iter<'a, D>(VecDeque<(&'a String, IterLapper<'a, u64, D>)>);

impl<'a, D> Iterator for Iter<'a, D> {
    type Item = (GenomicRange, &'a D);

    fn next(&mut self) -> Option<Self::Item> {
        if self.0.is_empty() {
            return None
        }
        match self.0[0].1.next() {
            None => {
                self.0.pop_front();
                self.next()
            }
            Some(interval) => {
                let bed = GenomicRange::new(self.0[0].0, interval.start, interval.stop);
                Some((bed, &interval.val))
            }
        }
    }
}

/// A Bed map that preserves the order.
#[derive(Debug, Clone)]
pub struct GIntervalIndexMap<D> {
    data: Vec<D>,
    indices: GIntervalMap<usize>,
}

impl<D> GIntervalIndexMap<D> {
    pub fn new() -> Self {
        Self { data: Vec::new(), indices: GIntervalMap::new() }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Note the results may be returned in arbitrary order.
    pub fn find<B: BEDLike>(&self, bed: &B) -> impl Iterator<Item = (GenomicRange, &D)> {
        self.indices.find(bed).map(|(g, i)| (g, &self.data[*i]))
    }

    /// Note the results may be returned in arbitrary order.
    pub fn find_index_of<B: BEDLike>(&self, bed: &B) -> GIntervalQueryIter<'_, usize> {
        self.indices.find(bed)
    }

    pub fn get(&self, index: usize) -> Option<&D> {
        self.data.get(index)
    }
}

impl<B: BEDLike, D> FromIterator<(B, D)> for GIntervalIndexMap<D> {
    fn from_iter<I: IntoIterator<Item = (B, D)>>(iter: I) -> Self {
        let mut data = Vec::new();
        let indices = iter.into_iter().enumerate().map(|(i, (bed, val))| {
            data.push(val);
            (bed, i)
        }).collect();
        Self { data, indices }
    }
}


/// A Bed map that preserves the order.
#[derive(Debug, Clone)]
pub struct GIntervalIndexSet {
    data: Vec<GenomicRange>,
    indices: GIntervalMap<usize>,
}

impl Index<usize> for GIntervalIndexSet {
    type Output = GenomicRange;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl GIntervalIndexSet {
    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn iter(&self) -> impl Iterator<Item = &GenomicRange> {
        self.data.iter()
    }

    pub fn is_overlapped<B: BEDLike>(&self, bed: &B) -> bool {
        self.find_full(bed).next().is_some()
    }

    /// Note the results may be returned in arbitrary order.
    pub fn find<B: BEDLike>(&self, bed: &B) -> impl Iterator<Item = GenomicRange> + '_ {
        self.indices.find(bed).map(|x| x.0)
    }

    /// Note the results may be returned in arbitrary order.
    pub fn find_index_of<B: BEDLike>(&self, bed: &B) -> impl Iterator<Item = usize> + '_ {
        self.indices.find(bed).map(|x| *x.1)
    }

    pub fn find_full<B: BEDLike>(&self, bed: &B) -> GIntervalQueryIter<'_, usize> {
        self.indices.find(bed)
    }

    pub fn get(&self, index: usize) -> Option<&GenomicRange> {
        self.data.get(index)
    }
}

impl IntoIterator for GIntervalIndexSet {
    type IntoIter = std::vec::IntoIter<Self::Item>;
    type Item = GenomicRange;

    fn into_iter(self) -> Self::IntoIter {
        self.data.into_iter()
    }
}

impl<B: BEDLike> FromIterator<B> for GIntervalIndexSet {
    fn from_iter<I: IntoIterator<Item = B>>(iter: I) -> Self {
        let mut data = Vec::new();
        let indices = iter.into_iter().enumerate().map(|(i, bed)| {
            data.push(bed.to_genomic_range());
            (bed, i)
        }).collect();
        Self { data, indices }
    }
}

#[cfg(test)]
mod bed_intersect_tests {
    use super::*;
    use itertools::Itertools;

    #[test]
    fn test_create() {
        let beds: Vec<GenomicRange> = [
            "chr1:200-500",
            "chr1:200-500",
            "chr1:1000-2000",
            "chr1:1000-2000",
        ].into_iter().map(|x| x.parse().unwrap()).collect();
        let values = vec![2, 3, 5, 5];
        let mut map: GIntervalMap<_> = beds.into_iter().zip(values.into_iter()).collect();
        map.insert(&GenomicRange::new("chr1", 1000, 2000), 5);
        map.insert(&GenomicRange::new("chr1", 1000, 2000), 6);

        let result: Vec<_> = map.iter().sorted().collect();
        let expected = vec![
            (GenomicRange::new("chr1", 200, 500), &2),
            (GenomicRange::new("chr1", 200, 500), &3),
            (GenomicRange::new("chr1", 1000, 2000), &5), 
            (GenomicRange::new("chr1", 1000, 2000), &5), 
            (GenomicRange::new("chr1", 1000, 2000), &5),
            (GenomicRange::new("chr1", 1000, 2000), &6),
        ];
        assert_eq!(result, expected);
    }

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
        let tree: GIntervalIndexSet = bed_set1.clone().into_iter().collect();

        assert_eq!(
            tree.find(&bed_set2[0]).collect::<Vec<_>>(),
            vec![bed_set1.as_slice()[0].clone(), bed_set1.as_slice()[0].clone()]
        );

        let result: Vec<GenomicRange> = bed_set2.into_iter().filter(|x| tree.is_overlapped(x)).collect();
        let expected = vec!["chr1:100-210".parse().unwrap()];
        assert_eq!(result, expected);

    }
}