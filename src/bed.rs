pub mod io;
pub mod map;

mod bed_trait;
pub use bed_trait::*;
mod score;
use bincode::{Decode, Encode};
pub use score::Score;
mod strand;
pub use strand::Strand;

use std::{fmt::{self, Write}, ops::Deref, str::FromStr};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

const DELIMITER: char = '\t';
const MISSING_ITEM : &str = ".";

/// A minimal BED record with only 3 fields.
#[derive(Encode, Decode, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GenomicRange(String, u64, u64);

impl GenomicRange {
    pub fn new<C>(chrom: C, start: u64, end: u64) -> Self
    where
        C: Into<String>,
    { Self(chrom.into(), start, end) }

    /// Convert the record to a string representation: chr:start-end
    pub fn pretty_show(&self) -> String {
        format!("{}:{}-{}", self.0, self.1, self.2)
    }
}

/// Convert string to GenomicRange. '\t', ':', and '-' are all considered as
/// valid delimiters. So any of the following formats is valid:
/// * chr1\t100\t200
/// * chr1:100-200
impl FromStr for GenomicRange {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut fields = s.split(&['\t', ':', '-']);
        let chrom = parse_chrom(&mut fields)?;
        let start = parse_start(&mut fields)?;
        let end = parse_end(&mut fields)?;
        Ok(GenomicRange::new(chrom, start, end))
    }
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

impl fmt::Display for GenomicRange {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}{}{}{}", self.chrom(), DELIMITER, self.start(),
            DELIMITER, self.end()
        )?;
        Ok(())
    }
}


/// A standard BED record.
#[derive(Encode, Decode, Clone, Debug, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
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

/// Generic BED record optional fields.
#[derive(Encode, Decode, Clone, Debug, Default, Eq, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
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


/// A NarrowPeak record is a BED6+4 format that is used to store called peaks.
#[derive(Encode, Decode, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct NarrowPeak {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: Option<String>,
    pub score: Option<Score>,
    pub strand: Option<Strand>,
    pub signal_value: f64,
    pub p_value: Option<f64>, 
    pub q_value: Option<f64>, 
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
            DELIMITER, self.p_value.unwrap_or(-1.0),
            DELIMITER, self.q_value.unwrap_or(-1.0),
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
            p_value: parse_pvalue(&mut fields).unwrap(),
            q_value: parse_pvalue(&mut fields).unwrap(),
            peak: fields.next().unwrap().parse().unwrap(),
        })
    }
}

/// A BroadPeak record is a BED6+4 format that is used to store called peaks.
#[derive(Encode, Decode, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BroadPeak {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub name: Option<String>,
    pub score: Option<Score>,
    pub strand: Option<Strand>,
    pub signal_value: f64,
    pub p_value: Option<f64>, 
    pub q_value: Option<f64>, 
}

impl BEDLike for BroadPeak {
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

impl fmt::Display for BroadPeak {
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
            "{}{}{}{}{}{}",
            DELIMITER, self.signal_value,
            DELIMITER, self.p_value.unwrap_or(-1.0),
            DELIMITER, self.q_value.unwrap_or(-1.0),
        )?;

        Ok(())
    }
}

impl FromStr for BroadPeak {
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
            p_value: parse_pvalue(&mut fields).unwrap(),
            q_value: parse_pvalue(&mut fields).unwrap(),
        })
    }
}

/// The bedGraph format allows display of continuous-valued data in track format.
/// This display type is useful for probability scores and transcriptome data. 
#[derive(Encode, Decode, Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BedGraph<V> {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub value: V,
}

impl<V> BedGraph<V> {
    pub fn new<C>(chrom: C, start: u64, end: u64, value: V) -> Self
    where
        C: Into<String>,
    { Self { chrom: chrom.into(), start, end, value } }

    pub fn from_bed<B: BEDLike>(bed: &B, value: V) -> Self {
        Self::new(bed.chrom(), bed.start(), bed.end(), value)
    }
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

fn parse_pvalue<'a, I>(fields: &mut I) -> Result<Option<f64>, ParseError>
where
    I: Iterator<Item = &'a str>,
{
    fields
        .next()
        .ok_or(ParseError::MissingScore)
        .and_then(|s| {
            let p = s.parse().unwrap();
            if p < 0.0 { Ok(None) } else { Ok(Some(p)) }
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

#[cfg(test)]
mod bed_tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let fields = OptionalFields::default();
        assert_eq!(fields.to_string(), "");

        let fields = OptionalFields::from(vec![String::from("n")]);
        assert_eq!(fields.to_string(), "n");

        let fields = OptionalFields::from(vec![String::from("n"), String::from("d")]);
        assert_eq!(fields.to_string(), "n\td");

        let genomic_range = GenomicRange::new("chr1", 100, 200);
        assert_eq!(genomic_range, GenomicRange::from_str("chr1\t100\t200").unwrap());
        assert_eq!(genomic_range, GenomicRange::from_str("chr1-100-200").unwrap());
        assert_eq!(genomic_range, GenomicRange::from_str("chr1:100-200").unwrap());
        assert_eq!(genomic_range, GenomicRange::from_str("chr1:100:200").unwrap());
        assert_eq!(genomic_range, GenomicRange::from_str(&genomic_range.pretty_show()).unwrap());
    }
}