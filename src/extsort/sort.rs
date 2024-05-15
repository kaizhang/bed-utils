use crate::extsort::merger::BinaryHeapMerger;
use crate::extsort::chunk::{ExternalChunk, ExternalChunkError};

use rayon::prelude::*;
use serde::de::DeserializeOwned;
use serde::Serialize;
use std::{fmt::{self, Display}, io};
use rayon;
use std::error::Error;
use bincode;
use std::{
    cmp::Ordering,
    path::{Path, PathBuf},
};

/// Sorting error.
#[derive(Debug)]
pub enum SortError<I: Error> {
    /// Temporary directory or file creation error.
    TempDir(io::Error),
    /// Workers thread pool initialization error.
    ThreadPoolBuildError(rayon::ThreadPoolBuildError),
    /// Common I/O error.
    IO(io::Error),
    /// Data serialization error.
    SerializationError(bincode::Error),
    /// Data deserialization error.
    DeserializationError(bincode::Error),
    /// Input data stream error
    InputError(I),
}

impl<I> Error for SortError<I>
where
    I: Error + 'static,
{
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(match &self {
            SortError::TempDir(err) => err,
            SortError::ThreadPoolBuildError(err) => err,
            SortError::IO(err) => err,
            SortError::SerializationError(err) => err,
            SortError::DeserializationError(err) => err,
            SortError::InputError(err) => err,
        })
    }
}

impl<I: Error> Display for SortError<I> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self {
            SortError::TempDir(err) => write!(f, "temporary directory or file not created: {}", err),
            SortError::ThreadPoolBuildError(err) => write!(f, "thread pool initialization failed: {}", err),
            SortError::IO(err) => write!(f, "I/O operation failed: {}", err),
            SortError::SerializationError(err) => write!(f, "data serialization error: {}", err),
            SortError::DeserializationError(err) => write!(f, "data deserialization error: {}", err),
            SortError::InputError(err) => write!(f, "input data stream error: {}", err),
        }
    }
}

/// Exposes external sorting (i.e. on disk sorting) capability on arbitrarily
/// sized iterator, even if the generated content of the iterator doesn't fit in
/// memory.
pub struct ExternalSorterBuilder {
    chunk_size: usize,
    tmp_dir: Option<PathBuf>,
    num_threads: Option<usize>,
}

impl ExternalSorterBuilder {
    pub fn new() -> Self {
        Self {
            chunk_size: 50000000,
            tmp_dir: None,
            num_threads: None,
        }
    }

    /// Sets the maximum size of each segment in number of sorted items.
    ///
    /// This number of items needs to fit in memory. While sorting, a
    /// in-memory buffer is used to collect the items to be sorted. Once
    /// it reaches the maximum size, it is sorted and then written to disk.
    ///
    /// Using a higher segment size makes sorting faster by leveraging
    /// faster in-memory operations.
    pub fn with_chunk_size(mut self, size: usize) -> Self {
        self.chunk_size = size;
        self
    }

    /// Sets directory in which sorted segments will be written (if it doesn't
    /// fit in memory).
    pub fn with_tmp_dir<P: AsRef<Path>>(mut self, path: P) -> Self {
        self.tmp_dir = Some(path.as_ref().to_path_buf());
        self
    }

    /// Uses Rayon to sort the in-memory buffer.
    ///
    /// This may not be needed if the buffer isn't big enough for parallelism to
    /// be gainful over the overhead of multithreading.
    pub fn num_threads(mut self, num_threads: usize) -> Self {
        self.num_threads = Some(num_threads);
        self
    }

    pub fn build(self) -> io::Result<ExternalSorter> {
        Ok(ExternalSorter {
            chunk_size: self.chunk_size,
            tmp_dir: _init_tmp_directory(self.tmp_dir.as_deref())?,
            thread_pool: _init_thread_pool(self.num_threads)?,
        })
    }
}

pub struct ExternalSorter {
    chunk_size: usize,
    /// Sorting thread pool.
    thread_pool: rayon::ThreadPool,
    /// Directory to be used to store temporary data.
    tmp_dir: tempfile::TempDir,
}

impl ExternalSorter {
    /// Sorts data from the input.
    /// Returns an iterator that can be used to get sorted data stream.
    ///
    /// # Arguments
    /// * `input` - Input stream data to be fetched from
    pub fn sort<I, T, E>(&self, input: I) -> Result<
        BinaryHeapMerger<T, ExternalChunkError, impl Fn(&T, &T) -> Ordering + Copy, ExternalChunk<T>>, SortError<E>
    >
    where
        T: Serialize + DeserializeOwned + Send + Ord,
        I: IntoIterator<Item = Result<T, E>>,
        E: Error,
    {
        self.sort_by(input, T::cmp)
    }

    /// Sorts a given iterator with a comparator function, returning a new iterator with items
    pub fn sort_by<I, T, E, F>(&self, input: I, cmp: F) -> Result<BinaryHeapMerger<T, ExternalChunkError, F, ExternalChunk<T>>, SortError<E>>
    where
        T: Serialize + DeserializeOwned + Send,
        I: IntoIterator<Item = Result<T, E>>,
        E: Error,
        F: Fn(&T, &T) -> Ordering + Sync + Send + Copy,
    {
        let mut chunk_buf = Vec::with_capacity(self.chunk_size);
        let mut external_chunks = Vec::new();

        for item in input.into_iter() {
            match item {
                Ok(item) => chunk_buf.push(item),
                Err(err) => return Err(SortError::InputError(err)),
            }

            if chunk_buf.len() >= self.chunk_size {
                external_chunks.push(self.create_chunk(chunk_buf, cmp)?);
                chunk_buf = Vec::with_capacity(self.chunk_size);
            }
        }

        if chunk_buf.len() > 0 {
            external_chunks.push(self.create_chunk(chunk_buf, cmp)?);
        }

        return Ok(BinaryHeapMerger::new(external_chunks, cmp));
    }

    fn create_chunk<T, F, E>(&self, mut buffer: Vec<T>, compare: F) -> Result<ExternalChunk<T>, SortError<E>>
    where
        T: Serialize + DeserializeOwned + Send,
        E: Error,
        F: Fn(&T, &T) -> Ordering + Sync + Send,
    {
        self.thread_pool.install(|| {
            buffer.par_sort_by(compare);
        });

        let external_chunk =
            ExternalChunk::build(&self.tmp_dir, buffer, None).map_err(|err| match err {
                ExternalChunkError::IO(err) => SortError::IO(err),
                ExternalChunkError::SerializationError(err) => SortError::SerializationError(err),
            })?;

        return Ok(external_chunk);
    }
}

fn _init_tmp_directory(
    tmp_path: Option<&Path>,
) -> io::Result<tempfile::TempDir> {
    if let Some(tmp_path) = tmp_path {
        tempfile::tempdir_in(tmp_path)
    } else {
        tempfile::tempdir()
    }
}

fn _init_thread_pool(
    threads_number: Option<usize>,
) -> io::Result<rayon::ThreadPool> {
    let mut thread_pool_builder = rayon::ThreadPoolBuilder::new();
    if let Some(threads_number) = threads_number {
        thread_pool_builder = thread_pool_builder.num_threads(threads_number);
    }
    thread_pool_builder.build().map_err(|x| io::Error::new(io::ErrorKind::Other, x))
}

#[cfg(test)]
mod test {
    use std::io;
    use std::path::Path;

    use rand::seq::SliceRandom;
    use rstest::*;

    use super::{ExternalSorter, ExternalSorterBuilder};

    #[rstest]
    #[case(false)]
    #[case(true)]
    fn test_external_sorter(#[case] reversed: bool) {
        let input_sorted = 0..100;

        let mut input: Vec<Result<i32, io::Error>> = Vec::from_iter(input_sorted.clone().map(|item| Ok(item)));
        input.shuffle(&mut rand::thread_rng());

        let sorter: ExternalSorter = ExternalSorterBuilder::new()
            .num_threads(2)
            .with_tmp_dir(Path::new("./"))
            .build()
            .unwrap();

        let compare = if reversed {
            |a: &i32, b: &i32| a.cmp(b).reverse()
        } else {
            |a: &i32, b: &i32| a.cmp(b)
        };

        let result = sorter.sort_by(input, compare).unwrap();

        let actual_result: Result<Vec<i32>, _> = result.collect();
        let actual_result = actual_result.unwrap();
        let expected_result = if reversed {
            Vec::from_iter(input_sorted.clone().rev())
        } else {
            Vec::from_iter(input_sorted.clone())
        };

        assert_eq!(actual_result, expected_result)
    }
}