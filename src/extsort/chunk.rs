//! External chunk.

use bitcode::{DecodeOwned, Encode};
use byteorder::{ReadBytesExt, WriteBytesExt};
use lz4::{Decoder, Encoder, EncoderBuilder};
use std::{
    error::Error,
    fmt::{self, Display},
    fs::File,
    io::{self, Read, Seek, Write},
    marker::PhantomData,
};

/// External chunk error
#[derive(Debug)]
pub enum ExternalChunkError {
    /// Common I/O error.
    IO(io::Error),
    /// Data serialization error.
    EncodeError(bitcode::Error),
}

impl From<io::Error> for ExternalChunkError {
    fn from(err: io::Error) -> Self {
        ExternalChunkError::IO(err)
    }
}

impl From<bitcode::Error> for ExternalChunkError {
    fn from(err: bitcode::Error) -> Self {
        ExternalChunkError::EncodeError(err)
    }
}

impl Error for ExternalChunkError {}

impl Display for ExternalChunkError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExternalChunkError::IO(err) => write!(f, "{}", err),
            ExternalChunkError::EncodeError(err) => write!(f, "{}", err),
        }
    }
}

/// External chunk interface. Provides methods for creating a chunk stored on file system and reading data from it.
pub struct ExternalChunk<T> {
    reader: Decoder<File>,
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
        file: File,
        items: impl IntoIterator<Item = T>,
        compression: u32,
    ) -> Result<Self, ExternalChunkError> {
        let mut builder = ExternalChunkBuilder::new(file, compression)?;
        for item in items.into_iter() {
            builder.add(item)?;
        }
        builder.finish()
    }
}

impl<T> Iterator for ExternalChunk<T>
where
    T: DecodeOwned,
{
    type Item = Result<T, ExternalChunkError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.read_u64::<byteorder::LittleEndian>() {
            Err(err) => match err.kind() {
                std::io::ErrorKind::UnexpectedEof => None,
                _ => Some(Err(ExternalChunkError::IO(err))),
            },
            Ok(length) => {
                let mut buf = vec![0u8; length as usize];
                if let Err(err) = self.reader.read_exact(buf.as_mut()) {
                    return Some(Err(ExternalChunkError::IO(err)));
                } else {
                    match bitcode::decode::<T>(&buf) {
                        Err(err) => Some(Err(ExternalChunkError::from(err))),
                        Ok(ser) => Some(Ok(ser)),
                    }
                }
            }
        }
    }
}

pub struct ExternalChunkBuilder<T> {
    writer: Encoder<File>,
    item_type: PhantomData<T>,
}

impl<T: Encode> ExternalChunkBuilder<T> {
    pub fn new(file: File, compression: u32) -> Result<Self, ExternalChunkError> {
        let writer = EncoderBuilder::new().level(compression).build(file)?;
        Ok(Self {
            writer,
            item_type: PhantomData,
        })
    }

    pub fn add(&mut self, item: T) -> Result<(), ExternalChunkError> {
        let result = bitcode::encode(&item);
        self.writer.write_u64::<byteorder::LittleEndian>(result.len() as u64)?;
        self.writer.write(&result)?;
        Ok(())
    }

    pub fn finish(self) -> Result<ExternalChunk<T>, ExternalChunkError> {
        let mut file = self.writer.finish().0;
        file.rewind()?;
        let reader = Decoder::new(file)?;

        Ok(ExternalChunk {
            reader,
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

        let tmp_file = tempfile::tempfile_in(&tmp_dir).unwrap();
        let chunk: ExternalChunk<i32> = ExternalChunk::new(tmp_file, saved.clone(), 0).unwrap();
        let restored = chunk.collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(restored, saved);

        let tmp_file = tempfile::tempfile_in(&tmp_dir).unwrap();
        let chunk: ExternalChunk<i32> =
            ExternalChunk::new(tmp_file, saved.clone(), 3).unwrap();
        let restored = chunk.collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(restored, saved);
    }
}
