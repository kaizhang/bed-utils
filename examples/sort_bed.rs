use std::{env, fs::File, io};
use bed_utils::bed::{io::Reader, BEDLike, BED};
use bed_utils::extsort::ExternalSorterBuilder;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");
    let file = File::open(src)?;
    let mut input = Reader::new(file, None);
    ExternalSorterBuilder::new()
        .with_chunk_size(5000000)
        .with_tmp_dir("./")
        .build().unwrap()
        .sort_by(input.records::<BED<5>>(), BEDLike::compare).unwrap()
        .for_each(|x| println!("{}", x.unwrap()));

    Ok(())
}