#![allow(dead_code)]

use urng_qmc::{QmcWord, QuasiRng};

/// SplitMix64, a fixed-seed generator for reproducible tests.
pub struct SplitMix64(pub u64);

impl SplitMix64 {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, bound)`.
    pub fn below(&mut self, bound: u64) -> u64 {
        ((self.next_u64() as u128 * bound as u128) >> 64) as u64
    }

    pub fn word<W: QmcWord>(&mut self) -> W {
        W::from_index(self.next_u64() >> (64 - W::BITS))
    }
}

impl urng_core::Rng for SplitMix64 {
    type Word = u64;
    fn nextu(&mut self) -> u64 {
        self.next_u64()
    }
}

/// `splitmix64_finalizer(h ^ x)`, the checksum step of tools/gen_vectors.py.
pub fn mix(h: u64, x: u64) -> u64 {
    let mut z = h ^ x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// The fixed-point word as a float in `[0, 1)`, computed exactly from the top 53 bits.
///
/// Independent of `urng_core::Word::to_float`, so the quality tests do not depend on it.
pub fn unit<W: QmcWord>(x: W) -> f64 {
    let top = x.to_u64() << (64 - W::BITS);
    (top >> 11) as f64 * (-53f64).exp2()
}

pub fn units<W: QmcWord, const D: usize>(p: [W; D]) -> [f64; D] {
    p.map(unit)
}

/// L2-star discrepancy by Warnock's formula, `O(N^2 D)`.
pub fn l2_star<const D: usize>(points: &[[f64; D]]) -> f64 {
    let n = points.len() as f64;
    let t1 = 3f64.powi(-(D as i32));
    let t2: f64 = points
        .iter()
        .map(|p| p.iter().map(|x| (1.0 - x * x) / 2.0).product::<f64>())
        .sum::<f64>()
        * 2.0
        / n;
    let mut t3 = 0.0;
    for p in points {
        for q in points {
            t3 += p
                .iter()
                .zip(q)
                .map(|(a, b)| 1.0 - a.max(*b))
                .product::<f64>();
        }
    }
    (t1 - t2 + t3 / (n * n)).sqrt()
}

/// Root of the expected squared L2-star discrepancy of `n` uniform random points.
pub fn l2_star_random(d: usize, n: usize) -> f64 {
    ((2f64.powi(-(d as i32)) - 3f64.powi(-(d as i32))) / n as f64).sqrt()
}

/// Least-squares slope of `ys` against `xs`.
pub fn slope(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let cov: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let var: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    cov / var
}

/// Top 16 bits of each coordinate of the first `2^m_max` points, one `Vec` per dimension.
pub fn top_bits<Q: QuasiRng<D>, const D: usize>(mut q: Q, m_max: u32) -> Vec<Vec<u16>> {
    let shift = Q::Word::BITS - 16;
    let mut dims = vec![Vec::with_capacity(1 << m_max); D];
    for _ in 0..1u32 << m_max {
        for (d, x) in q.next_point().into_iter().enumerate() {
            dims[d].push((x.to_u64() >> shift) as u16);
        }
    }
    dims
}

/// For every `m <= m_max`, the first `2^m` values have distinct top `m` bits: one value in
/// each interval `[k / 2^m, (k + 1) / 2^m)`.
pub fn assert_1d_stratified(values: &[u16], m_max: u32, what: &str) {
    for m in 1..=m_max {
        let mut seen = vec![false; 1 << m];
        for &v in &values[..1 << m] {
            let cell = (v >> (16 - m)) as usize;
            assert!(!seen[cell], "{what}: m={m}, two points in cell {cell}");
            seen[cell] = true;
        }
    }
}

/// Dimensions 1 and 2 form a (0, 2)-sequence: every block of `2^m` points starting at a
/// multiple of `2^m` has one point in each `2^a x 2^(m-a)` cell, for every `a`.
pub fn assert_2d_t0(x: &[u16], y: &[u16], m_max: u32, what: &str) {
    for m in 1..=m_max {
        let blocks = (x.len() >> m).min(4);
        for block in 0..blocks {
            let range = block << m..(block + 1) << m;
            for a in 0..=m {
                let mut seen = vec![false; 1 << m];
                for i in range.clone() {
                    let cx = (u32::from(x[i]) >> (16 - a)) as usize;
                    let cy = (u32::from(y[i]) >> (16 - (m - a))) as usize;
                    let cell = (cx << (m - a)) | cy;
                    assert!(!seen[cell], "{what}: m={m} a={a} block={block}");
                    seen[cell] = true;
                }
            }
        }
    }
}
