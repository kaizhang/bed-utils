use std::{env, fs::File, io};
use extsort::sorter::ExternalSorter;
use bed_utils::bed::{BEDLike, BED, io::Reader};
use std::path::PathBuf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");
    let file = File::open(src)?;
    let sorter = ExternalSorter::new()
        .with_segment_size(50000000)
        .with_sort_dir(PathBuf::from("./"))
        .with_parallel_sort();
    sorter.sort_by(
        Reader::new(file, None).records::<BED<6>>().map(|x| x.unwrap()),
        BEDLike::compare,
    )?.for_each(|x| println!("{}", x));

    Ok(())
}