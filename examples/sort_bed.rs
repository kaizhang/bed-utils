use std::{env, fs::File, io};
use bed_utils::bed::{io::Reader, sort_bed, BEDLike, BED};
use std::path::PathBuf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");
    let file = File::open(src)?;
    sort_bed(
        Reader::new(file, None).records::<BED<5>>(),
        Some("./"),
    ).for_each(|x| println!("{}", x.unwrap()));

    Ok(())
}