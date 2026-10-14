mod common;

use urng_core::Rng;
use urng_qmc::{
    AsRng, DigitalNet, DigitalShifted, QuasiRng, Rd32, Rd64, Seek, Shifted, Sobol32, Sobol64,
    VanDerCorput64, Weyl64,
};

use common::{SplitMix64, assert_2d_t0, top_bits};

/// Returns 1, 2, 3, ...
struct Counter32(u32);

impl Rng for Counter32 {
    type Word = u32;
    fn nextu(&mut self) -> u32 {
        self.0 += 1;
        self.0
    }
}

#[test]
fn shifted_adds_the_shift() {
    let mut rng = SplitMix64(1);
    let mut q = Shifted::new(Rd64::<3>::new(), &mut rng);
    let shift = *q.shift();
    let mut plain = Rd64::<3>::new();
    for _ in 0..1000 {
        let p = plain.next_point();
        assert_eq!(
            q.next_point(),
            core::array::from_fn(|j| p[j].wrapping_add(shift[j]))
        );
    }
    assert_eq!(q.get_ref().index(), 1000);
}

#[test]
fn digital_shifted_xors_the_mask() {
    let mut rng = SplitMix64(2);
    let mut q = DigitalShifted::new(Sobol32::<5>::new(), &mut rng);
    let mask = *q.mask();
    let mut plain = Sobol32::<5>::new();
    for _ in 0..1000 {
        let p = plain.next_point();
        assert_eq!(q.next_point(), core::array::from_fn(|j| p[j] ^ mask[j]));
    }
    assert_eq!(q.into_inner(), plain);
}

#[test]
fn shift_words_come_from_the_rng() {
    // 64-bit sequence, 32-bit generator: two draws per word, the first one high.
    let q = DigitalShifted::new(Sobol64::<2>::new(), &mut Counter32(0));
    assert_eq!(q.mask(), &[0x0000_0001_0000_0002, 0x0000_0003_0000_0004]);
    // 32-bit sequence, 64-bit generator: the high half of each draw.
    let mut rng = SplitMix64(3);
    let mut copy = SplitMix64(3);
    let q = Shifted::new(Rd32::<2>::new(), &mut rng);
    let want: [u32; 2] = [
        (copy.next_u64() >> 32) as u32,
        (copy.next_u64() >> 32) as u32,
    ];
    assert_eq!(q.shift(), &want);
}

#[test]
fn zero_shift_is_the_identity() {
    let mut a = Shifted::with_shift(Weyl64::golden(), [0]);
    let mut b = DigitalShifted::with_mask(VanDerCorput64::new(), [0]);
    let mut weyl = Weyl64::golden();
    let mut vdc = VanDerCorput64::new();
    for _ in 0..100 {
        assert_eq!(a.next_point(), weyl.next_point());
        assert_eq!(b.next_point(), vdc.next_point());
    }
}

#[test]
fn fill_point_matches_next_point() {
    let mut rng = SplitMix64(4);
    let a = DigitalShifted::new(Sobol64::<16>::new(), &mut rng);
    let (mut a, mut b) = (a.clone(), a);
    let mut out = [0; 16];
    for _ in 0..100 {
        b.fill_point(&mut out);
        assert_eq!(a.next_point(), out);
    }
}

/// XOR-ing a digital net keeps every elementary interval filled exactly once.
#[test]
fn digital_shift_keeps_the_net_property() {
    for seed in 0..10 {
        let mut rng = SplitMix64(seed);
        let d64 = top_bits(DigitalShifted::new(Sobol64::<2>::new(), &mut rng), 14);
        assert_2d_t0(&d64[0], &d64[1], 12, &format!("Sobol64 seed {seed}"));
        let d32 = top_bits(DigitalShifted::new(Sobol32::<2>::new(), &mut rng), 14);
        assert_2d_t0(&d32[0], &d32[1], 12, &format!("Sobol32 seed {seed}"));
    }
}

/// A digitally shifted net is still a digital net, so shifts compose.
#[test]
fn digital_shifts_compose() {
    fn is_net<Q: DigitalNet<D>, const D: usize>(_: &Q) {}
    let inner = DigitalShifted::with_mask(Sobol64::<2>::new(), [1, 2]);
    let mut q = DigitalShifted::with_mask(inner, [4, 8]);
    is_net(&q);
    let mut plain = Sobol64::<2>::new();
    let _ = (q.next_point(), plain.next_point());
    let p = plain.next_point();
    assert_eq!(q.next_point(), [p[0] ^ 5, p[1] ^ 10]);
}

#[test]
fn one_dimensional_wrapper_works_as_rng() {
    let mut rng = AsRng::new(Shifted::with_shift(Weyl64::golden(), [1]));
    assert_eq!(rng.nextu(), 1);
    assert_eq!(rng.nextu(), 1u64.wrapping_add(Weyl64::golden().alpha()));
    rng.seek(0);
    assert_eq!(rng.nextu(), 1);
}

#[test]
fn every_coordinate_is_uniform_after_shifting() {
    // The first point of an unscrambled sequence is the origin; after a random shift its
    // coordinates follow the shift, so the mean over many shifts is close to 1/2.
    let mut rng = SplitMix64(5);
    let n = 4096;
    let mean = (0..n)
        .map(|_| {
            let mut q = DigitalShifted::new(Sobol64::<1>::new(), &mut rng);
            common::unit(q.next_point()[0])
        })
        .sum::<f64>()
        / n as f64;
    assert!((mean - 0.5).abs() < 0.02, "mean {mean}");
}
