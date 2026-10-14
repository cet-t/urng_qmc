#![allow(dead_code)]

use urng_qmc::QmcWord;

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

/// The fixed-point word as a float in `[0, 1)`, computed exactly from the top 53 bits.
///
/// The tests convert with this rather than `urng_core::Word::to_float`, whose `u64`
/// version drops the top bit in urng_core 0.1.1.
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
