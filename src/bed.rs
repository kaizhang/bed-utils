pub mod io;
pub mod tree;

pub use self::{score::Score, strand::Strand};
use itertools::Itertools;
mod score;
mod strand;

use std::{
    fmt::{self, Write},
    ops::Deref,
    str::FromStr,
    cmp::Ordering,
};

const DELIMITER: char = '\t';
const MISSING_ITEM : &str = ".";

/// A minimal BED record with only 3 fields.
#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct GenomicRange(String, u64, u64);

impl GenomicRange {
    pub fn new<C>(chrom: C, start: u64, end: u64) -> Self
    where
        C: Into<String>,
    { Self(chrom.into(), start, end) }
}

impl FromStr for GenomicRange {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let chrom = parse_chrom(&mut fields)?;
        let start = parse_start(&mut fields)?;
        let end = parse_end(&mut fields)?;
        Ok(GenomicRange::new(chrom, start, end))
    }
}

/// A standard BED record.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BED<const N: u8> {
    chrom: String,
    start: u64,
    end: u64,
    pub name: Option<String>,
    pub score: Option<Score>,
    pub strand: Option<Strand>,
    pub optional_fields: OptionalFields,
}

impl<const N: u8> BED<N> {
    pub fn new<C>(chrom: C, start: u64, end: u64, name: Option<String>,
        score: Option<Score>, strand: Option<Strand>, optional_fields: OptionalFields) -> Self
    where
        C: Into<String>,
    { Self { chrom: chrom.into(), start, end, name, score, strand, optional_fields } }
}

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

    /// Convert the record to a string representation: chr:start-end
    fn to_string(&self) -> String {
        format!("{}:{}-{}", self.chrom(), self.start(), self.end())
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

pub fn merge_bed_with<I, B, O, F>(bed_iter: I, merger: F) -> MergeBed<std::vec::IntoIter<B>, B, F>
where
    I: Iterator<Item = B>,
    B: BEDLike,
    F: Fn(Vec<B>) -> O,
{
    merge_sorted_bed_with(bed_iter.sorted_by(BEDLike::compare), merger)
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
) -> MergeBed<std::vec::IntoIter<B>, B, impl Fn(Vec<B>) -> GenomicRange>
where
    I: Iterator<Item = B>,
    B: BEDLike,
{
    merge_sorted_bed(bed_iter.sorted_by(BEDLike::compare))
}



impl BEDLike for GenomicRange {
    fn chrom(&self) -> &str { &self.0 }
    fn set_chrom(&mut self, chrom: &str) -> &mut Self {
        self.0 = chrom.to_string();
        self
    }
    fn start(&self) -> u64 { self.1 }
    fn set_start(&mut self, start: u64) -> &mut Self {
        self.1 = start;
        self
    }
    fn end(&self) -> u64 { self.2 }
    fn set_end(&mut self, end: u64) -> &mut Self {
        self.2 = end;
        self
    }
    fn name(&self) -> Option<&str> { None }
    fn score(&self) -> Option<Score> { None }
    fn strand(&self) -> Option<Strand> { None }
}

impl<const N: u8> BEDLike for BED<N> {
    fn chrom(&self) -> &str { &self.chrom }
    fn set_chrom(&mut self, chrom: &str) -> &mut Self {
        self.chrom = chrom.to_string();
        self
    }
    fn start(&self) -> u64 { self.start }
    fn set_start(&mut self, start: u64) -> &mut Self {
        self.start = start;
        self
    }
    fn end(&self) -> u64 { self.end }
    fn set_end(&mut self, end: u64) -> &mut Self {
        self.end = end;
        self
    }
    fn name(&self) -> Option<&str> { self.name.as_deref() }
    fn score(&self) -> Option<Score> { self.score }
    fn strand(&self) -> Option<Strand> { self.strand }
}

// Display trait
impl<const N: u8> fmt::Display for BED<N> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{}{}{}{}",
            self.chrom(),
            DELIMITER,
            self.start(),
            DELIMITER,
            self.end()
        )?;
        if N > 3 {
            write!(f, "{}{}", DELIMITER, self.name().unwrap_or(MISSING_ITEM))?;
            if N > 4 {
                f.write_char(DELIMITER)?;
                if let Some(score) = self.score() {
                    write!(f, "{}", score)?;
                } else { f.write_str(MISSING_ITEM)?; }

                if N > 5 {
                    f.write_char(DELIMITER)?;
                    if let Some(strand) = self.strand() {
                        write!(f, "{}", strand)?;
                    } else { f.write_str(MISSING_ITEM)?; }
                }
            }
        }
        Ok(())
    }
}

impl<const N: u8> FromStr for BED<N> {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        let chrom = parse_chrom(&mut fields)?;
        let start = parse_start(&mut fields)?;
        let end = parse_end(&mut fields)?;
        let name = if N > 3 { parse_name(&mut fields)? } else { None };
        let score = if N > 4 { parse_score(&mut fields)? } else { None };
        let strand = if N > 5 { parse_strand(&mut fields)? } else { None };
        Ok(BED::new(chrom, start, end, name, score, strand, OptionalFields::default()))
    }
}

fn parse_chrom<'a, I>(fields: &mut I) -> Result<&'a str, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingReferenceSequenceName)
}

fn parse_start<'a, I>(fields: &mut I) -> Result<u64, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingStartPosition)
        .and_then(|s| lexical::parse(s).map_err(ParseError::InvalidStartPosition))
}

fn parse_end<'a, I>(fields: &mut I) -> Result<u64, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingEndPosition)
        .and_then(|s| lexical::parse(s).map_err(ParseError::InvalidEndPosition))
}

fn parse_name<'a, I>(fields: &mut I) -> Result<Option<String>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingName)
        .map(|s| match s {
            MISSING_ITEM => None,
            _ => Some(s.into()),
        })
}

fn parse_score<'a, I>(fields: &mut I) -> Result<Option<Score>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingScore)
        .and_then(|s| match s {
            MISSING_ITEM => Ok(None),
            _ => s.parse().map(Some).map_err(ParseError::InvalidScore),
        })
}

fn parse_strand<'a, I>(fields: &mut I) -> Result<Option<Strand>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingStrand)
        .and_then(|s| match s {
            MISSING_ITEM => Ok(None),
            _ => s.parse().map(Some).map_err(ParseError::InvalidStrand),
        })
}

/// An error returned when a raw BED record fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The reference sequence name is missing.
    MissingReferenceSequenceName,
    /// The start position is missing.
    MissingStartPosition,
    /// The start position is invalid.
    InvalidStartPosition(lexical::Error),
    /// The end position is missing.
    MissingEndPosition,
    /// The end position is invalid.
    InvalidEndPosition(lexical::Error),
    /// The name is missing.
    MissingName,
    /// The score is missing.
    MissingScore,
    /// The score is invalid.
    InvalidScore(score::ParseError),
    /// The strand is missing.
    MissingStrand,
    /// The strand is invalid.
    InvalidStrand(strand::ParseError),
}


/// Generic BED record optional fields.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct OptionalFields(Vec<String>);

impl Deref for OptionalFields {
    type Target = [String];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl fmt::Display for OptionalFields {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, field) in self.0.iter().enumerate() {
            if i > 0 {
                f.write_char(DELIMITER)?;
            }

            f.write_str(field)?;
        }

        Ok(())
    }
}

impl From<Vec<String>> for OptionalFields {
    fn from(fields: Vec<String>) -> Self {
        Self(fields)
    }
}


/// A standard BED record.
#[derive(Clone, Debug, PartialEq)]
pub struct NarrowPeak {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: Option<String>,
    pub score: Option<Score>,
    pub strand: Option<Strand>,
    pub signal_value: f64,
    pub p_value: f64, 
    pub q_value: f64, 
    pub peak: u64, 
}

impl BEDLike for NarrowPeak {
    fn chrom(&self) -> &str { &self.chrom }
    fn set_chrom(&mut self, chrom: &str) -> &mut Self {
        self.chrom = chrom.to_string();
        self
    }
    fn start(&self) -> u64 { self.start }
    fn set_start(&mut self, start: u64) -> &mut Self {
        self.start = start;
        self
    }
    fn end(&self) -> u64 { self.end }
    fn set_end(&mut self, end: u64) -> &mut Self {
        self.end = end;
        self
    }
    fn name(&self) -> Option<&str> { self.name.as_deref() }
    fn score(&self) -> Option<Score> { self.score }
    fn strand(&self) -> Option<Strand> { self.strand }
}

impl fmt::Display for NarrowPeak {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{}{}{}{}{}{}",
            self.chrom(),
            DELIMITER, self.start(),
            DELIMITER, self.end(),
            DELIMITER, self.name().unwrap_or(MISSING_ITEM),
        )?;

        f.write_char(DELIMITER)?;
        if let Some(x) = self.score() {
            write!(f, "{}", x)?;
        } else {
            f.write_str(MISSING_ITEM)?;
        }
        f.write_char(DELIMITER)?;
        if let Some(x) = self.strand() {
            write!(f, "{}", x)?;
        } else {
            f.write_str(MISSING_ITEM)?;
        }
        write!(
            f,
            "{}{}{}{}{}{}{}{}",
            DELIMITER, self.signal_value,
            DELIMITER, self.p_value,
            DELIMITER, self.q_value,
            DELIMITER, self.peak,
        )?;

        Ok(())
    }
}

impl FromStr for NarrowPeak {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(DELIMITER);
        Ok(Self {
            chrom: parse_chrom(&mut fields)?.to_string(),
            start: parse_start(&mut fields)?,
            end: parse_end(&mut fields)?,
            name: parse_name(&mut fields)?,
            score: parse_score(&mut fields)?,
            strand: parse_strand(&mut fields)?,
            signal_value: fields.next().unwrap().parse().unwrap(),
            p_value: fields.next().unwrap().parse().unwrap(),
            q_value: fields.next().unwrap().parse().unwrap(),
            peak: fields.next().unwrap().parse().unwrap(),
        })
    }
}

/// The bedGraph format allows display of continuous-valued data in track format.
/// This display type is useful for probability scores and transcriptome data. 
#[derive(Clone, Debug, PartialEq)]
pub struct BedGraph<V> {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub value: V,
}

impl<V> BEDLike for BedGraph<V> {
    fn chrom(&self) -> &str { &self.chrom }
    fn set_chrom(&mut self, chrom: &str) -> &mut Self {
        self.chrom = chrom.to_string();
        self
    }
    fn start(&self) -> u64 { self.start }
    fn set_start(&mut self, start: u64) -> &mut Self {
        self.start = start;
        self
    }
    fn end(&self) -> u64 { self.end }
    fn set_end(&mut self, end: u64) -> &mut Self {
        self.end = end;
        self
    }
    fn name(&self) -> Option<&str> { None }
    fn score(&self) -> Option<Score> { None }
    fn strand(&self) -> Option<Strand> { None }
}

impl<V> fmt::Display for BedGraph<V>
where
    V: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result
    {
        write!(
            f,
            "{}{}{}{}{}{}{}",
            self.chrom(),
            DELIMITER, self.start(),
            DELIMITER, self.end(),
            DELIMITER, self.value,
        )
    }
}

impl<V> FromStr for BedGraph<V>
where
    V: FromStr,
    <V as FromStr>::Err: std::fmt::Debug,
{
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err>
    {
        let mut fields = s.split(DELIMITER);
        Ok(Self {
            chrom: parse_chrom(&mut fields)?.to_string(),
            start: parse_start(&mut fields)?,
            end: parse_end(&mut fields)?,
            value: fields.next().unwrap().parse().unwrap(),
        })
    }
}


#[cfg(test)]
mod bed_tests {
    use super::*;

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
            (500, 1000)
        ].into_iter().map(|(a,b)| GenomicRange::new("chr1", a, b));
        let expect: Vec<GenomicRange> = [
            (0, 220),
            (500, 1000)
        ].into_iter().map(|(a,b)| GenomicRange::new("chr1", a, b)).collect();
        assert_eq!(merge_sorted_bed(input.clone()).collect::<Vec<_>>(), expect);
        assert_eq!(merge_bed(input.rev()).collect::<Vec<_>>(), expect);
    }

    #[test]
    fn test_fmt() {
        let fields = OptionalFields::default();
        assert_eq!(fields.to_string(), "");

        let fields = OptionalFields::from(vec![String::from("n")]);
        assert_eq!(fields.to_string(), "n");

        let fields = OptionalFields::from(vec![String::from("n"), String::from("d")]);
        assert_eq!(fields.to_string(), "n\td");
    }
}