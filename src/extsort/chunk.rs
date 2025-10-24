//! External chunk.

use bincode::{Decode, Encode};
use lz4::{Decoder, Encoder, EncoderBuilder};
use std::error::Error;
use std::fmt::{self, Display};
use std::fs::File;
use std::io::{self, BufWriter};
use std::io::{prelude::*, BufReader};
use std::marker::PhantomData;

use byteorder::{ReadBytesExt, WriteBytesExt};
use tempfile;

/// External chunk error
#[derive(Debug)]
pub enum ExternalChunkError {
    /// Common I/O error.
    IO(io::Error),
    /// Data serialization error.
    EncodeError(bincode::error::EncodeError),
    DecodeError(bincode::error::DecodeError),
}

impl From<io::Error> for ExternalChunkError {
    fn from(err: io::Error) -> Self {
        ExternalChunkError::IO(err)
    }
}

impl From<bincode::error::EncodeError> for ExternalChunkError {
    fn from(err: bincode::error::EncodeError) -> Self {
        ExternalChunkError::EncodeError(err)
    }
}

impl From<bincode::error::DecodeError> for ExternalChunkError {
    fn from(err: bincode::error::DecodeError) -> Self {
        ExternalChunkError::DecodeError(err)
    }
}

impl Error for ExternalChunkError {}

impl Display for ExternalChunkError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExternalChunkError::IO(err) => write!(f, "{}", err),
            ExternalChunkError::EncodeError(err) => write!(f, "{}", err),
            ExternalChunkError::DecodeError(err) => write!(f, "{}", err),
        }
    }
}

/// External chunk interface. Provides methods for creating a chunk stored on file system and reading data from it.
pub struct ExternalChunk<T> {
    reader: Box<dyn Read>,
    item_type: PhantomData<T>,
}

impl<T> ExternalChunk<T>
where
    T: Encode,
{
    /// Builds an instance of an external chunk creating file and dumping the items to it.
    ///
    /// # Arguments
    /// * `dir` - Directory the chunk file is created in
    /// * `items` - Items to be dumped to the chunk
    /// * `buf_size` - File I/O buffer size
    pub(crate) fn new(
        dir: &tempfile::TempDir,
        items: impl IntoIterator<Item = T>,
        compression: Option<u32>,
    ) -> Result<Self, ExternalChunkError> {
        let mut builder = ExternalChunkBuilder::new(dir, compression)?;
        for item in items.into_iter() {
            builder.add(item)?;
        }
        builder.finish()
    }
}

impl<T> Iterator for ExternalChunk<T>
where
    T: Decode<()>,
{
    type Item = Result<T, ExternalChunkError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_u64::<byteorder::LittleEndian>() {
            Err(err) => match err.kind() {
                std::io::ErrorKind::UnexpectedEof => None,
                _ => Some(Err(ExternalChunkError::IO(err))),
            },
            Ok(length) => {
                let config = bincode::config::standard();
                let mut buf = vec![0u8; length as usize];
                if let Err(err) = self.reader.read_exact(buf.as_mut()) {
                    return Some(Err(ExternalChunkError::IO(err)));
                } else {
                    match bincode::decode_from_slice(&buf, config) {
                        Err(err) => Some(Err(ExternalChunkError::from(err))),
                        Ok((ser, n)) => {
                            if n != length as usize {
                                Some(Err(ExternalChunkError::IO(io::Error::new(
                                    io::ErrorKind::InvalidData,
                                    format!("Expected {} bytes, got {}", length, n),
                                ))))
                            } else {
                                Some(Ok(ser))
                            }
                        }
                    }
                }
            }
        }
    }
}

pub struct ExternalChunkBuilder<T> {
    writer: Result<Encoder<File>, BufWriter<File>>,
    item_type: PhantomData<T>,
}

impl<T: Encode> ExternalChunkBuilder<T> {
    pub fn new(
        dir: &tempfile::TempDir,
        compression: Option<u32>,
    ) -> Result<Self, ExternalChunkError> {
        let tmp_file = tempfile::tempfile_in(dir)?;

        let writer = if let Some(lvl) = compression {
            Ok(EncoderBuilder::new().level(lvl).build(tmp_file)?)
        } else {
            Err(BufWriter::new(tmp_file))
        };

        Ok(Self {
            writer,
            item_type: PhantomData,
        })
    }

    pub fn add(&mut self, item: T) -> Result<(), ExternalChunkError> {
        let result = bincode::encode_to_vec(&item, bincode::config::standard())
            .map_err(ExternalChunkError::from)?;

        self.writer.as_mut().map_or_else(
            |w| {
                w.write_u64::<byteorder::LittleEndian>(result.len() as u64)?;
                w.write(&result)?;
                Ok(())
            },
            |w| {
                w.write_u64::<byteorder::LittleEndian>(result.len() as u64)?;
                w.write(&result)?;
                Ok(())
            },
        )
    }

    pub fn finish(self) -> Result<ExternalChunk<T>, ExternalChunkError> {
        let reader: Result<Box<dyn Read>, ExternalChunkError> = self.writer.map_or_else(
            |w| {
                let mut file = w.into_inner().unwrap();
                file.rewind()?;
                let reader: Box<dyn Read> = Box::new(BufReader::new(file));
                Ok(reader)
            },
            |w| {
                let mut file = w.finish().0;
                file.rewind()?;
                let reader = Box::new(Decoder::new(file)?);
                Ok(reader)
            },
        );

        Ok(ExternalChunk {
            reader: reader?,
            item_type: PhantomData,
        })
    }
}

#[cfg(test)]
mod test {
    use rstest::*;

    use super::ExternalChunk;

    #[fixture]
    fn tmp_dir() -> tempfile::TempDir {
        tempfile::tempdir_in("./").unwrap()
    }

    #[rstest]
    fn test_chunk(tmp_dir: tempfile::TempDir) {
        let saved = Vec::from_iter(0..100);

        let chunk: ExternalChunk<i32> = ExternalChunk::new(&tmp_dir, saved.clone(), None).unwrap();
        let restored = chunk.collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(restored, saved);

        let chunk: ExternalChunk<i32> =
            ExternalChunk::new(&tmp_dir, saved.clone(), Some(3)).unwrap();
        let restored = chunk.collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(restored, saved);
    }
}
