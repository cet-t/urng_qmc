mod common;
#[allow(dead_code)]
#[path = "data/sobol_vectors.rs"]
mod vectors;

use urng_qmc::{
    QuasiRng, RD_MAX_DIM, Rd32, Rd64, SOBOL_MAX_DIM, Sobol32, Sobol64, VanDerCorput32,
    VanDerCorput64, Weyl32, Weyl64,
};

use common::unit;

#[test]
fn weyl_golden() {
    assert_eq!(Weyl64::golden().alpha(), 0x9E37_79B9_7F4A_7C15);
    assert_eq!(Weyl32::golden().alpha(), 0x9E37_79B9);
    assert_eq!(Weyl64::default(), Weyl64::golden());

    let mut q = Weyl64::golden();
    for n in 0u64..100 {
        assert_eq!(q.next_point(), [n.wrapping_mul(0x9E37_79B9_7F4A_7C15)]);
    }
}

#[test]
fn weyl_makes_alpha_odd() {
    assert_eq!(Weyl32::new(4).alpha(), 5);
    assert_eq!(Weyl64::new(0).alpha(), 1);
    assert_eq!(Weyl64::new(7).alpha(), 7);
}

#[test]
fn weyl_seed_shifts_every_point() {
    let mut plain = Weyl32::golden();
    let mut seeded = Weyl32::golden().with_seed(0x8000_0000);
    assert_eq!(seeded.seed(), 0x8000_0000);
    for _ in 0..100 {
        let [x] = plain.next_point();
        assert_eq!(seeded.next_point(), [x.wrapping_add(0x8000_0000)]);
    }
    // Seeding later keeps the index.
    let late = plain.clone().with_seed(3);
    let [x] = plain.next_point();
    assert_eq!(late.clone().next_point(), [x.wrapping_add(3)]);
}

#[test]
fn rd_constants_match_spec() {
    assert_eq!(Rd64::<1>::ALPHA, [0x9E37_79B9_7F4A_7C15]);
    assert_eq!(
        Rd64::<2>::ALPHA,
        [0xC13F_A9A9_02A6_328F, 0x91E1_0DA5_C79E_7B1D]
    );
    assert_eq!(
        Rd64::<3>::ALPHA,
        [
            0xD1B5_4A32_D192_ED03,
            0xABC9_8388_FB8F_AC03,
            0x8CB9_2BA7_2F3D_8DD7
        ]
    );
    assert_eq!(Rd32::<1>::ALPHA, [0x9E37_79B9]);
    assert_eq!(Rd32::<2>::ALPHA, [0xC13F_A9A9, 0x91E1_0DA5]);
    assert_eq!(Rd32::<3>::ALPHA, [0xD1B5_4A33, 0xABC9_8389, 0x8CB9_2BA7]);
    assert_eq!(Weyl64::golden().alpha(), Rd64::<1>::ALPHA[0]);
    assert_eq!(Weyl32::golden().alpha(), Rd32::<1>::ALPHA[0]);
}

/// Positive root of `x^(d+1) = x + 1` in f64.
fn phi(d: i32) -> f64 {
    let mut x = 2.0f64;
    for _ in 0..100 {
        x -= (x.powi(d + 1) - x - 1.0) / ((d + 1) as f64 * x.powi(d) - 1.0);
    }
    x
}

fn check_rd<const D: usize>() {
    let p = phi(D as i32);
    for (j, (&a64, &a32)) in Rd64::<D>::ALPHA.iter().zip(&Rd32::<D>::ALPHA).enumerate() {
        assert_eq!(a64 & 1, 1, "D={D} j={j}");
        assert_eq!(a32 & 1, 1, "D={D} j={j}");
        assert_eq!(a32, (a64 >> 32) as u32 | 1, "D={D} j={j}");
        // f64 carries 53 bits and powi rounds once per factor, so this only catches gross
        // errors such as a wrong row; tools/gen_rd.py checks the exact values.
        let approx = p.powi(-(j as i32 + 1)) * 2f64.powi(64);
        let rel = (a64 as f64 - approx).abs() / approx;
        assert!(rel <= 2f64.powi(-48), "D={D} j={j}: relative error {rel:e}");
    }
}

macro_rules! check_rd_dims {
    ($($d:literal)+) => { $(check_rd::<$d>();)+ };
}

#[test]
fn rd_constants_are_odd_and_close_to_theory() {
    assert_eq!(RD_MAX_DIM, 16);
    check_rd_dims!(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16);
}

#[test]
fn rd_first_points() {
    let mut q = Rd64::<3>::new();
    for n in 0u64..100 {
        assert_eq!(q.next_point(), Rd64::<3>::ALPHA.map(|a| a.wrapping_mul(n)));
    }
    let mut q = Rd32::<2>::new().with_seed([1, 2]);
    for n in 0u32..100 {
        assert_eq!(
            q.next_point(),
            [
                1u32.wrapping_add(Rd32::<2>::ALPHA[0].wrapping_mul(n)),
                2u32.wrapping_add(Rd32::<2>::ALPHA[1].wrapping_mul(n)),
            ]
        );
    }
}

#[test]
fn van_der_corput_first_points() {
    let expected = [0.0, 0.5, 0.25, 0.75, 0.125, 0.625, 0.375, 0.875];
    let mut q32 = VanDerCorput32::new();
    let mut q64 = VanDerCorput64::new();
    for want in expected {
        assert_eq!(unit(q32.next_point()[0]), want);
        assert_eq!(unit(q64.next_point()[0]), want);
    }
    let mut q = VanDerCorput64::new();
    for n in 0u64..1000 {
        assert_eq!(q.next_point(), [n.reverse_bits()]);
    }
}

#[test]
fn sobol_spec_reference_points() {
    let expected = [
        (0.0, 0.0),
        (0.5, 0.5),
        (0.75, 0.25),
        (0.25, 0.75),
        (0.375, 0.375),
        (0.875, 0.875),
        (0.625, 0.125),
        (0.125, 0.625),
    ];
    let mut q32 = Sobol32::<2>::new();
    let mut q64 = Sobol64::<2>::new();
    for (x, y) in expected {
        let [a, b] = q32.next_point();
        assert_eq!((unit(a), unit(b)), (x, y));
        let [a, b] = q64.next_point();
        assert_eq!((unit(a), unit(b)), (x, y));
    }
    // f32 conversion through urng_core is exact for these points.
    let mut q = Sobol32::<2>::new();
    for (x, y) in expected {
        assert_eq!(q.next_point_f(), [x as f32, y as f32]);
    }
}

#[test]
#[ignore = "urng_core 0.1.1 drops the top bit in u64::to_float; re-enable once urng_core 0.1.2 is required"]
fn sobol64_next_point_f() {
    let mut q = Sobol64::<2>::new();
    let _ = q.next_point();
    assert_eq!(q.next_point_f(), [0.5, 0.5]);
    assert_eq!(q.next_point_f(), [0.75, 0.25]);
}

#[test]
fn sobol_matches_scipy_d8() {
    let mut q32 = Sobol32::<8>::new();
    let mut q64 = Sobol64::<8>::new();
    for n in 0..vectors::POINTS {
        assert_eq!(
            q32.next_point(),
            vectors::SOBOL32_D8[n],
            "Sobol32 point {n}"
        );
        assert_eq!(
            q64.next_point(),
            vectors::SOBOL64_D8[n],
            "Sobol64 point {n}"
        );
    }
}

#[test]
fn sobol_matches_scipy_d1024_checksums() {
    assert_eq!(SOBOL_MAX_DIM, 1024);
    let mut q = Sobol32::<1024>::new();
    let mut sums = [0u32; 1024];
    let mut p = [0u32; 1024];
    for _ in 0..vectors::POINTS {
        q.fill_point(&mut p);
        for (h, x) in sums.iter_mut().zip(p) {
            *h = h.rotate_left(5) ^ x;
        }
    }
    for (d, (got, want)) in sums
        .iter()
        .zip(&vectors::SOBOL32_D1024_CHECKSUM)
        .enumerate()
    {
        assert_eq!(got, want, "dimension {}", d + 1);
    }
}

#[test]
fn sobol64_high_half_is_sobol32() {
    let mut q32 = Sobol32::<1024>::new();
    let mut q64 = Sobol64::<1024>::new();
    let (mut p32, mut p64) = ([0u32; 1024], [0u64; 1024]);
    for n in 0..4096 {
        q32.fill_point(&mut p32);
        q64.fill_point(&mut p64);
        assert_eq!(p64.map(|x| (x >> 32) as u32), p32, "point {n}");
    }
}

#[test]
fn fill_point_matches_next_point() {
    let mut a = Sobol64::<5>::new();
    let mut b = Sobol64::<5>::new();
    let mut p = [0; 5];
    for _ in 0..100 {
        b.fill_point(&mut p);
        assert_eq!(a.next_point(), p);
    }
    assert_eq!(a, b);
}
