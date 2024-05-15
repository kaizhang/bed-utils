//! External chunk.

use std::error::Error;
use std::fmt::{self, Display};
use std::io::{self, BufReader, BufWriter};
use std::io::prelude::*;
use std::marker::PhantomData;
use lz4::{Decoder, EncoderBuilder};

use bincode::Options;
use byteorder::{ReadBytesExt, WriteBytesExt};
use serde::Serialize;
use tempfile;

/// External chunk error
#[derive(Debug)]
pub enum ExternalChunkError {
    /// Common I/O error.
    IO(io::Error),
    /// Data serialization error.
    SerializationError(bincode::Error),
}

impl From<io::Error> for ExternalChunkError {
    fn from(err: io::Error) -> Self {
        ExternalChunkError::IO(err)
    }
}

impl From<bincode::Error> for ExternalChunkError {
    fn from(err: bincode::Error) -> Self {
        ExternalChunkError::SerializationError(err)
    }
}

impl Error for ExternalChunkError {}

impl Display for ExternalChunkError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExternalChunkError::IO(err) => write!(f, "{}", err),
            ExternalChunkError::SerializationError(err) => write!(f, "{}", err),
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
    T: serde::ser::Serialize + serde::de::DeserializeOwned,
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
        let mut tmp_file = tempfile::tempfile_in(dir)?;

        if let Some(level) = compression {
            let mut writer = EncoderBuilder::new().level(level).build(tmp_file.try_clone()?)?;
            dump(&mut writer, items)?;
            writer.finish().1?;
        } else {
            let mut writer = BufWriter::new(tmp_file.try_clone()?);
            dump(&mut writer, items)?;
            writer.flush()?;
        }

        tmp_file.rewind()?;
        if let Some(_) = compression {
            let reader = Box::new(Decoder::new(tmp_file)?);
            Ok(Self { reader, item_type: PhantomData })
        } else {
            let reader = Box::new(BufReader::new(tmp_file));
            Ok(Self { reader, item_type: PhantomData })
        }
    }
}

fn dump<T: Serialize, R: Write>(
    chunk_writer: &mut R,
    items: impl IntoIterator<Item = T>,
) -> Result<(), ExternalChunkError> {
    let ser = bincode::DefaultOptions::new();
    for item in items.into_iter() {
        let result = ser.serialize(&item).map_err(ExternalChunkError::from)?;
        let size = result.len() as u64;
        chunk_writer.write_u64::<byteorder::LittleEndian>(size)?;
        chunk_writer.write(result.as_slice())?;
    }
    return Ok(());
}

impl<T> Iterator for ExternalChunk<T>
where
    T: serde::ser::Serialize + serde::de::DeserializeOwned,
{
    type Item = Result<T, ExternalChunkError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_u64::<byteorder::LittleEndian>() {
            Err(err) => match err.kind() {
                std::io::ErrorKind::UnexpectedEof => None,
                _ => Some(Err(ExternalChunkError::IO(err))),
            },
            Ok(size) => {
                let mut buf = vec![0u8; size as usize];
                if let Err(err) =  self.reader.read_exact(buf.as_mut()) {
                    return Some(Err(ExternalChunkError::IO(err)));
                } else {
                    let ser = bincode::DefaultOptions::new();
                    Some(ser.deserialize(&buf).map_err(ExternalChunkError::from))
                }
            }
        }
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

        let restored: Result<Vec<i32>, _> = chunk.collect();
        let restored = restored.unwrap();

        assert_eq!(restored, saved);
    }
}