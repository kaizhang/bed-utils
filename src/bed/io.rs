use crate::bed::BEDLike;

use std::{
    io::{self, Error, ErrorKind, Read, Write, BufRead, BufReader},
};
use std::str::FromStr;
use std::marker::PhantomData;

use super::ParseError;

/// An iterator over records of a FASTQ reader.
///
/// This is created by calling [`Reader::records`].
pub struct Records<'a, B, R> {
    inner: &'a mut Reader<R>,
    buf: String,
    phantom: PhantomData<B>,
}

impl<'a, B, R> Records<'a, B, R>
where
    R: Read,
    B: FromStr,
{
    pub fn new(inner: &'a mut Reader<R>) -> Self {
        Self {
            inner,
            buf: String::new(),
            phantom: PhantomData,
        }
    }
}

impl<'a, B, R> Iterator for Records<'a, B, R>
where
    R: Read,
    B: FromStr<Err = ParseError>,
{
    type Item = io::Result<B>;

    fn next(&mut self) -> Option<Self::Item> {
        self.buf.clear();
        match self.inner.read_record(&mut self.buf) {
            Ok(LineSize::Size(0)) => None,
            Ok(LineSize::Skip) => self.next(),
            Ok(_) => Some(self.buf.parse().map_err(
                |e| Error::new(ErrorKind::Other, format!("{:?}: {}", e, &self.buf))
            )),
            Err(e) => Some(Err(e)),
        }
    }
}

/// An iterator over records of a FASTQ reader.
///
/// This is created by calling [`Reader::records`].
pub struct IntoRecords<B, R> {
    inner: Reader<R>,
    buf: String,
    phantom: PhantomData<B>,
}

impl<B, R> IntoRecords<B, R>
where
    R: Read,
    B: FromStr,
{
    pub fn new(inner: Reader<R>) -> Self {
        Self { inner, buf: String::new(), phantom: PhantomData }
    }
}

impl<B, R> Iterator for IntoRecords<B, R>
where
    R: Read,
    B: FromStr<Err = ParseError>,
{
    type Item = io::Result<B>;

    fn next(&mut self) -> Option<Self::Item> {
        self.buf.clear();
        match self.inner.read_record(&mut self.buf) {
            Ok(LineSize::Size(0)) => None,
            Ok(LineSize::Skip) => self.next(),
            Ok(_) => Some(self.buf.parse().map_err(
                |e| Error::new(ErrorKind::Other, format!("{:?}: {}", e, &self.buf))
            )),
            Err(e) => Some(Err(e)),
        }
    }
}

/// A BED reader.
pub struct Reader<R> {
    inner: BufReader<R>,
    skip_start_with: Option<String>,
}

impl<R> Reader<R>
where
    R: Read,
{
    /// Creates a BED reader.
    pub fn new(inner: R, skip_start_with: Option<String>) -> Self {
        Self { inner: BufReader::new(inner), skip_start_with }
    }

    /// Reads a single raw BED record.
    pub fn read_record(&mut self, buf: &mut String) -> io::Result<LineSize> {
        let size = read_line(&mut self.inner, buf)?;
        if size > 0 && self.skip_start_with.as_ref().map_or(false, |x| buf.starts_with(x)) {
            Ok(LineSize::Skip)
        } else {
            Ok(LineSize::Size(size))
        }
    }

    /// Returns an iterator over records starting from the current stream position.
    ///
    /// The stream is expected to be at the start of a record.
    ///
    pub fn records<B: FromStr + BEDLike>(&mut self) -> Records<'_, B, R> {
        Records::new(self)
    }

    pub fn into_records<B: FromStr + BEDLike>(self) -> IntoRecords<B, R> {
        IntoRecords::new(self)
    }
}

pub enum LineSize {
    Size(usize),
    Skip,
}

fn read_line<R>(reader: &mut R, buf: &mut String) -> io::Result<usize>
where
    R: BufRead,
{
    const LINE_FEED: char = '\n';
    const CARRIAGE_RETURN: char = '\r';

    match reader.read_line(buf) {
        Ok(0) => Ok(0),
        Ok(n) => {
            if buf.ends_with(LINE_FEED) {
                buf.pop();

                if buf.ends_with(CARRIAGE_RETURN) {
                    buf.pop();
                }
            }
            Ok(n)
        }
        Err(e) => Err(e),
    }
}

/// A BED writer.
pub struct Writer<W> { inner: W }

impl<W> Writer<W> where W: Write {
    /// Creates a BED writer.
    pub fn new(inner: W) -> Self { Self { inner } }

    /// Writes a BED record.
    pub fn write_record<B>(&mut self, record: &B) -> io::Result<()>
    where
        B: std::fmt::Display + BEDLike,
    {
        writeln!(&mut self.inner, "{}", record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::BED;

    #[test]
    fn test_read_line() {
        fn t(buf: &mut String, mut reader: &[u8], expected: &str) {
            buf.clear();
            read_line(&mut reader, buf).unwrap();
            assert_eq!(buf, expected);
        }

        let mut buf = String::new();

        t(&mut buf, b"noodles\n", "noodles");
        t(&mut buf, b"noodles\r\n", "noodles");
        t(&mut buf, b"noodles", "noodles");
    }

    #[test]
    fn test_read_record() {
        let data = b"\
chr1	200	1000	r1	100	+
chr2	220	2000	r2	2	-
chr10	2000	10000	r3	3	+
" as &[u8];
        let mut reader = Reader::new(data, None);
        for b in reader.records() {
            let x: BED<6> = b.unwrap();
            println!("{}", x);
        }

        /*
        read_record(&mut reader, &mut record)?;
        assert_eq!(record, Record::new("noodles:1/1", "AGCT", "abcd"));

        read_record(&mut reader, &mut record)?;
        assert_eq!(record, Record::new("noodles:2/1", "TCGA", "dcba"));

        let n = read_record(&mut reader, &mut record)?;
        assert_eq!(n, 0);
        */

    }

    #[test]
    fn test_write_record() -> Result<(), Box<dyn std::error::Error>> {
        let mut writer = Writer::new(Vec::new());
        let record: BED<3> = "sq0\t8\t13".parse().unwrap();
        writer.write_record(&record)?;
        assert_eq!(writer.inner, b"sq0\t8\t13\n");
        Ok(())
    }
}


