//! Time per point for each sequence at D = 1, 2, 8, 32 (where supported).

use std::hint::black_box;
use std::time::Duration;

use criterion::{
    BenchmarkGroup, Criterion, criterion_group, criterion_main, measurement::WallTime,
};
use urng_qmc::{
    DigitalShifted, Halton64, QuasiRng, Rd64, Seek, Shifted, Sobol32, Sobol64, VanDerCorput64,
    Weyl64,
};

fn point<Q: QuasiRng<D>, const D: usize>(g: &mut BenchmarkGroup<WallTime>, name: &str, mut q: Q) {
    g.bench_function(format!("{name}/D={D}"), |b| {
        b.iter(|| black_box(q.next_point()))
    });
}

fn next_point(c: &mut Criterion) {
    let mut g = c.benchmark_group("next_point");
    point(&mut g, "Weyl64", Weyl64::golden());
    point(&mut g, "VanDerCorput64", VanDerCorput64::new());
    point(&mut g, "Rd64", Rd64::<1>::new());
    point(&mut g, "Rd64", Rd64::<2>::new());
    point(&mut g, "Rd64", Rd64::<8>::new());
    point(&mut g, "Sobol64", Sobol64::<1>::new());
    point(&mut g, "Sobol64", Sobol64::<2>::new());
    point(&mut g, "Sobol64", Sobol64::<8>::new());
    point(&mut g, "Sobol64", Sobol64::<32>::new());
    point(&mut g, "Sobol32", Sobol32::<32>::new());
    point(&mut g, "Halton64", Halton64::<1>::new());
    point(&mut g, "Halton64", Halton64::<2>::new());
    point(&mut g, "Halton64", Halton64::<8>::new());
    point(&mut g, "Halton64", Halton64::<32>::new());
    point(
        &mut g,
        "Shifted<Rd64>",
        Shifted::with_shift(Rd64::<8>::new(), [1; 8]),
    );
    point(
        &mut g,
        "DigitalShifted<Sobol64>",
        DigitalShifted::with_mask(Sobol64::<32>::new(), [1; 32]),
    );
    g.finish();
}

fn seek(c: &mut Criterion) {
    let mut g = c.benchmark_group("seek");
    let mut q = Sobol64::<32>::new();
    let mut n = 0x9E37_79B9_7F4A_7C15u64;
    g.bench_function("Sobol64/D=32", |b| {
        b.iter(|| {
            n = n.wrapping_add(0x9E37_79B9);
            q.seek(black_box(n));
        })
    });
    g.finish();
}

criterion_group! {
    name = benches;
    config = Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(2));
    targets = next_point, seek
}
criterion_main!(benches);
