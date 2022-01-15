mod perf;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fxd::FixedDecimal;

fn batch_from_str(ss: &[String], out: &mut [FixedDecimal]) {
    for (s, fd) in ss.iter().zip(out.iter_mut()) {
        fd.from_ascii_str(s, true).unwrap();
    }
}

fn bench_from_str(c: &mut Criterion) {
    use rand::prelude::*;
    let mut rng = rand::thread_rng();
    const N: usize = 1024;
    let ss: Vec<String> = (0..N)
        .into_iter()
        .map(|_| {
            let intg: u32 = rng.gen_range(0..1 << 30);
            let frac: u32 = rng.gen_range(0..1 << 20);
            if frac > 0 {
                format!("{}.{}", intg, frac)
            } else {
                format!("{}", intg)
            }
            // "Infinity".to_string()
        })
        .collect();

    let mut ds: Vec<FixedDecimal> = (0..N).into_iter().map(|_| FixedDecimal::zero()).collect();
    c.bench_function("batch-from-str", |b| {
        b.iter(|| batch_from_str(black_box(&ss), black_box(&mut ds)))
    });
}

criterion_group!(benches, bench_from_str);
criterion_main!(benches);
