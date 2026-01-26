//! BED record score.

use std::{error, fmt, str::FromStr};
use std::ops::Deref;

use bitcode::{Decode, Encode};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A BED record score.
#[derive(Encode, Decode, Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Score(u16);

impl fmt::Display for Score {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Deref for Score {
    type Target = u16;
    fn deref(&self) -> &Self::Target { &self.0 }
}

/// An error returned when a raw BED record score fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The input failed to be parsed as an integer.
    Parse(lexical::Error),
    /// The input is invalid.
    Invalid(TryFromIntError),
}

impl error::Error for ParseError {}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(e) => write!(f, "parse error: {}", e),
            Self::Invalid(e) => write!(f, "invalid input: {}", e),
        }
    }
}

impl FromStr for Score {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let n: u32 = lexical::parse(s).map_err(ParseError::Parse)?;
        Ok(Self::try_from(n).unwrap_or(Score(1000)))
        //map_err(ParseError::Invalid)
    }
}

/// An error returned when a raw BED record score fails to convert.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TryFromIntError(u32);

impl error::Error for TryFromIntError {}

impl fmt::Display for TryFromIntError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "invalid value: {}", self.0)
    }
}

impl TryFrom<u32> for Score {
    type Error = TryFromIntError;

    fn try_from(n: u32) -> Result<Self, Self::Error> {
        if n > 1000 {
            Err(TryFromIntError(n))
        } else {
            Ok(Self(n as u16))
        }
    }
}

impl TryFrom<u16> for Score {
    type Error = TryFromIntError;

    fn try_from(n: u16) -> Result<Self, Self::Error> {
        Self::try_from(n as u32)
    }
}

impl From<Score> for u16 {
    fn from(score: Score) -> Self {
        score.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        assert_eq!(Score(1).to_string(), "1");
        assert_eq!(Score(1000).to_string(), "1000");
    }

    #[test]
    fn test_try_from_u16_for_score() {
        assert_eq!(Score::try_from(1u16), Ok(Score(1)));
        assert_eq!(Score::try_from(1000u16), Ok(Score(1000)));
    }

    #[test]
    fn test_from_score_for_u16() {
        assert_eq!(u16::from(Score(8)), 8);
    }
}