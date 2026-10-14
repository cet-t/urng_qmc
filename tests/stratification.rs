use urng_qmc::{QmcWord, QuasiRng, Sobol32, Sobol64, VanDerCorput32, VanDerCorput64};

/// Top 16 bits of each coordinate of the first `2^m_max` points, one `Vec` per dimension.
fn top_bits<Q: QuasiRng<D>, const D: usize>(mut q: Q, m_max: u32) -> Vec<Vec<u16>> {
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
fn assert_1d_stratified(values: &[u16], m_max: u32, what: &str) {
    for m in 1..=m_max {
        let mut seen = vec![false; 1 << m];
        for &v in &values[..1 << m] {
            let cell = (v >> (16 - m)) as usize;
            assert!(!seen[cell], "{what}: m={m}, two points in cell {cell}");
            seen[cell] = true;
        }
    }
}

#[test]
fn van_der_corput_1d() {
    for (name, values) in [
        ("vdc32", top_bits(VanDerCorput32::new(), 16)),
        ("vdc64", top_bits(VanDerCorput64::new(), 16)),
    ] {
        assert_1d_stratified(&values[0], 16, name);
    }
}

#[test]
fn sobol_1d_first_64_dims() {
    for (bits, dims) in [
        (32, top_bits(Sobol32::<64>::new(), 16)),
        (64, top_bits(Sobol64::<64>::new(), 16)),
    ] {
        for (d, values) in dims.iter().enumerate() {
            assert_1d_stratified(values, 16, &format!("Sobol{bits} dimension {}", d + 1));
        }
    }
}

#[test]
fn sobol_1d_all_1024_dims() {
    for (d, values) in top_bits(Sobol64::<1024>::new(), 12).iter().enumerate() {
        assert_1d_stratified(values, 12, &format!("dimension {}", d + 1));
    }
}

/// Dimensions 1 and 2 form a (0, 2)-sequence: every block of `2^m` points starting at a
/// multiple of `2^m` has one point in each `2^a x 2^(m-a)` cell, for every `a`.
fn assert_2d_t0(x: &[u16], y: &[u16], m_max: u32, what: &str) {
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

#[test]
fn sobol_2d_t0() {
    let d32 = top_bits(Sobol32::<2>::new(), 14);
    assert_2d_t0(&d32[0], &d32[1], 12, "Sobol32");
    let d64 = top_bits(Sobol64::<2>::new(), 14);
    assert_2d_t0(&d64[0], &d64[1], 12, "Sobol64");
}
