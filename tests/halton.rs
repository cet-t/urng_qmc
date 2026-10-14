#[allow(dead_code)]
#[path = "data/halton_vectors.rs"]
mod vectors;

use urng_qmc::{
    HALTON_MAX_DIM, Halton32, Halton64, QuasiRng, Seek, VanDerCorput32, VanDerCorput64,
};

const PRIMES: [u64; 32] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131,
];

#[test]
fn matches_exact_values() {
    assert_eq!(HALTON_MAX_DIM, 32);
    let mut q = Halton64::<32>::new();
    for (n, want) in vectors::HALTON64_N.iter().zip(&vectors::HALTON64_D32) {
        q.seek(*n);
        assert_eq!(&q.next_point(), want, "Halton64 index {n}");
    }
    let mut q = Halton32::<32>::new();
    for (n, want) in vectors::HALTON32_N.iter().zip(&vectors::HALTON32_D32) {
        q.seek(*n);
        assert_eq!(&q.next_point(), want, "Halton32 index {n}");
    }
}

#[test]
fn first_coordinate_is_van_der_corput() {
    let (mut h32, mut h64) = (Halton32::<3>::new(), Halton64::<3>::new());
    let (mut v32, mut v64) = (VanDerCorput32::new(), VanDerCorput64::new());
    for _ in 0..10_000 {
        assert_eq!(h32.next_point()[0], v32.next_point()[0]);
        assert_eq!(h64.next_point()[0], v64.next_point()[0]);
    }
}

/// The first `b^m` points of coordinate `j` sit on the grid `i / b^m`, one per point.
fn check_grid<Q: QuasiRng<32>>(mut q: Q, bits: u32, word: impl Fn(Q::Word) -> u64) {
    let points: Vec<[u64; 32]> = (0..4096).map(|_| q.next_point().map(&word)).collect();
    for (j, &b) in PRIMES.iter().enumerate() {
        let mut bm = b;
        while bm <= 4096 {
            let mut seen = vec![false; bm as usize];
            for p in &points[..bm as usize] {
                let scaled = u128::from(p[j]) * u128::from(bm);
                // Nearest grid index, and how far the point is from it in output units.
                let i = (scaled + (1 << (bits - 1))) >> bits;
                let off = scaled.abs_diff(i << bits);
                assert!(off <= u128::from(bm), "base {b}: off grid by {off}");
                assert!(
                    !seen[i as usize],
                    "base {b}, b^m = {bm}: grid point {i} twice"
                );
                seen[i as usize] = true;
            }
            bm *= b;
        }
    }
}

#[test]
fn first_b_pow_m_points_fill_the_grid() {
    check_grid(Halton64::<32>::new(), 64, |x| x);
    check_grid(Halton32::<32>::new(), 32, u64::from);
}
