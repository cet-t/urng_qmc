//! `SobolWide32` / `SobolWide64` with all 21201 Joe–Kuo dimensions (`sobol-21201` feature).
#![cfg(feature = "sobol-21201")]

mod common;
#[allow(dead_code)]
#[path = "data/sobol_vectors.rs"]
mod vectors;

use urng_qmc::{QuasiRng, SOBOL_WIDE_MAX_DIM, Seek, Sobol32, Sobol64, SobolWide32, SobolWide64};

const D: usize = 21201;

/// Heap-allocated so the 170 KB points stay off the test thread's stack.
fn zeros<T: Copy + Default + std::fmt::Debug>() -> Box<[T; D]> {
    vec![T::default(); D].into_boxed_slice().try_into().unwrap()
}

#[test]
fn matches_scipy_checksum() {
    assert_eq!(SOBOL_WIDE_MAX_DIM, D);
    let mut q = Box::new(SobolWide32::<D>::new());
    let mut sums = zeros::<u64>();
    let mut p = zeros::<u32>();
    for _ in 0..vectors::POINTS {
        q.fill_point(&mut p);
        for (h, &x) in sums.iter_mut().zip(p.iter()) {
            *h = common::mix(*h, x.into());
        }
    }
    let folded = sums.iter().fold(0, |f, &h| common::mix(f, h));
    assert_eq!(folded, vectors::SOBOL32_D21201_FOLDED);
}

#[test]
fn high_half_is_sobol32_and_every_dimension_is_stratified() {
    let mut q32 = Box::new(SobolWide32::<D>::new());
    let mut q64 = Box::new(SobolWide64::<D>::new());
    let (mut p32, mut p64) = (zeros::<u32>(), zeros::<u64>());
    // Bit `x >> 24` of seen[d] is set once a point of dimension d fell in that cell.
    let mut seen = vec![[0u128; 2]; D];
    for n in 0..256 {
        q32.fill_point(&mut p32);
        q64.fill_point(&mut p64);
        for d in 0..D {
            assert_eq!(
                (p64[d] >> 32) as u32,
                p32[d],
                "point {n}, dimension {}",
                d + 1
            );
            let cell = (p32[d] >> 24) as usize;
            let bit = 1u128 << (cell % 128);
            assert_eq!(
                seen[d][cell / 128] & bit,
                0,
                "dimension {}: cell {cell} twice",
                d + 1
            );
            seen[d][cell / 128] |= bit;
        }
    }
}

#[test]
fn same_points_as_the_narrow_types() {
    let (mut w32, mut n32) = (SobolWide32::<1024>::new(), Sobol32::<1024>::new());
    let (mut w64, mut n64) = (SobolWide64::<1024>::new(), Sobol64::<1024>::new());
    for _ in 0..4096 {
        assert_eq!(w32.next_point(), n32.next_point());
        assert_eq!(w64.next_point(), n64.next_point());
    }
}

#[test]
fn seek_matches_sequential_across_chunks() {
    let mut walked = SobolWide64::<600>::new();
    for n in 0..5000u64 {
        let mut jumped = SobolWide64::<600>::new();
        jumped.seek(n);
        assert_eq!(jumped, walked, "seek({n})");
        let _ = walked.next_point();
    }
    let mut q = SobolWide32::<300>::new();
    q.seek(u64::from(u32::MAX));
    let _ = q.next_point();
    assert_eq!(q, SobolWide32::<300>::new());
}
