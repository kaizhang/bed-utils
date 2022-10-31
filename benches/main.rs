use bed_utils::bed;
use criterion::{criterion_group, criterion_main, Criterion};
use std::fs::File;
use std::io::Read;

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("IO");
    group.sample_size(30);

    let mut input = String::new();
    File::open("data/test.bed").unwrap().read_to_string(&mut input).unwrap();

    group.bench_with_input("IO", &input, |b, i| {
        b.iter(|| {
            bed::io::Reader::new(i.as_bytes(), None).into_records().for_each(|x| {
                let b: bed::BED<6> = x.unwrap();
                ()
            });
        });
    });

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);