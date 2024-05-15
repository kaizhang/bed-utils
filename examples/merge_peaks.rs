use std::{env, error::Error, fs::File, io};
use bed_utils::bed::{BEDLike, io::Reader, NarrowPeak, merge_bed_with};
use std::path::Path;
use flate2::read::MultiGzDecoder;

pub fn merge_peaks<I, E>(peaks: I, half_window_size: u64) -> impl Iterator<Item = Vec<NarrowPeak>>
where
    I: Iterator<Item = Result<NarrowPeak, E>>,
    E: Error,
{
    fn iterative_merge(mut peaks: Vec<NarrowPeak>) -> Vec<NarrowPeak> {
        let mut result = Vec::new();
        while !peaks.is_empty() {
            let best_peak = peaks.iter()
                .max_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap()).unwrap()
                .clone();
            peaks = peaks.into_iter().filter(|x| x.n_overlap(&best_peak) == 0).collect();
            result.push(best_peak);
        }
        result
    }

    merge_bed_with(
        peaks.map(move |r| r.map(|mut x| {
            let summit = x.start() + x.peak;
            x.start = summit.saturating_sub(half_window_size);
            x.end = summit + half_window_size + 1;
            x.peak = summit - x.start;
            x
        })),
        iterative_merge,
        None::<&str>,
    )
}

pub(crate) fn open_file<P: AsRef<Path>>(file: P) -> Box<dyn std::io::Read> {
    if is_gzipped(file.as_ref()) {
        Box::new(MultiGzDecoder::new(File::open(file.as_ref()).unwrap()))
    } else {
        Box::new(File::open(file.as_ref()).unwrap())
    }
}

/// Determine if a file is gzipped.
pub(crate) fn is_gzipped<P: AsRef<Path>>(file: P) -> bool {
    MultiGzDecoder::new(File::open(file).unwrap()).header().is_some()
}


fn main() -> io::Result<()> {
    let files: Vec<String> = env::args().skip(1).collect();

    let peak_iter = files.into_iter().flat_map(|fl|
        Reader::new(open_file(fl), None)
        .into_records::<NarrowPeak>()
    );
    let peaks: Vec<_> = merge_peaks(peak_iter, 250).flatten().collect();
    for peak in peaks {
        println!("{}", peak);
    }

    Ok(())
}