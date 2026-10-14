mod common;

use std::collections::BTreeSet;
use std::fmt::Debug;

use urng_qmc::{
    AsRng, QmcWord, QuasiRng, Rd32, Rd64, Seek, Sobol32, Sobol64, VanDerCorput32, VanDerCorput64,
    Weyl32, Weyl64,
};

use common::SplitMix64;

/// `seek(n)` gives the same state as `n` sequential updates, for every interesting `n`
/// below `limit`.
fn seek_matches_sequential<Q, const D: usize>(fresh: &Q, limit: u64)
where
    Q: QuasiRng<D> + Seek + Clone + PartialEq + Debug,
{
    let mut targets = BTreeSet::new();
    targets.extend([0, 1]);
    for k in 1..64 {
        for n in [(1u64 << k) - 1, 1 << k] {
            if n < limit {
                targets.insert(n);
            }
        }
    }
    let mut rng = SplitMix64(limit);
    for _ in 0..32 {
        targets.insert(rng.below(limit));
    }

    let mut walked = fresh.clone();
    let mut n = 0;
    for target in targets {
        while n < target {
            let _ = walked.next_point();
            n += 1;
        }
        let mut jumped = fresh.clone();
        jumped.seek(target);
        assert_eq!(jumped, walked, "seek({target})");
        assert_eq!(jumped.index(), target);
        assert_eq!(jumped.clone().next_point(), walked.clone().next_point());
    }
}

/// `seek(n)` then one update equals `seek(n + 1)`, for indices far beyond a walk.
fn seek_chain<Q, const D: usize>(fresh: &Q)
where
    Q: QuasiRng<D> + Seek + Clone + PartialEq + Debug,
{
    let bits = Q::Word::BITS;
    let mask = u64::MAX >> (64 - bits);
    let mut rng = SplitMix64(u64::from(bits) + D as u64);
    for _ in 0..100 {
        let n = rng.next_u64() & mask;
        let mut a = fresh.clone();
        a.seek(n);
        let _ = a.next_point();
        let mut b = fresh.clone();
        b.seek(n.wrapping_add(1) & mask);
        assert_eq!(a, b, "seek({n}) + 1");
    }
}

/// The last point is followed by point 0, and indices wrap modulo `2^BITS`.
fn seek_period<Q, const D: usize>(fresh: &Q)
where
    Q: QuasiRng<D> + Seek + Clone + PartialEq + Debug,
{
    let bits = Q::Word::BITS;
    let mask = u64::MAX >> (64 - bits);
    let mut q = fresh.clone();
    q.seek(mask);
    assert_eq!(q.index(), mask);
    let _ = q.next_point();
    assert_eq!(&q, fresh, "point 2^{bits} is point 0");
    assert_eq!(q.index(), 0);

    if bits < 64 {
        let mut a = fresh.clone();
        a.seek((1 << bits) + 5);
        let mut b = fresh.clone();
        b.seek(5);
        assert_eq!(a, b);
        assert_eq!(a.index(), 5);
    }
}

fn all<Q, const D: usize>(fresh: Q, limit: u64)
where
    Q: QuasiRng<D> + Seek + Clone + PartialEq + Debug,
{
    seek_matches_sequential(&fresh, limit);
    seek_chain(&fresh);
    seek_period(&fresh);
}

const LIMIT: u64 = 1 << 16;

#[test]
fn weyl() {
    all(Weyl32::golden(), LIMIT);
    all(Weyl64::golden(), LIMIT);
    all(Weyl32::new(12345).with_seed(0xDEAD_BEEF), LIMIT);
    all(Weyl64::golden().with_seed(u64::MAX), LIMIT);
}

#[test]
fn rd() {
    all(Rd32::<1>::new(), LIMIT);
    all(Rd32::<3>::new().with_seed([1, 2, 3]), LIMIT);
    all(Rd64::<2>::new(), LIMIT);
    all(Rd64::<16>::new().with_seed([u64::MAX; 16]), LIMIT);
}

#[test]
fn van_der_corput() {
    all(VanDerCorput32::new(), LIMIT);
    all(VanDerCorput64::new(), LIMIT);
}

#[test]
fn sobol() {
    all(Sobol32::<1>::new(), LIMIT);
    all(Sobol32::<8>::new(), LIMIT);
    all(Sobol64::<8>::new(), LIMIT);
    all(Sobol32::<1024>::new(), 1 << 10);
    all(Sobol64::<1024>::new(), 1 << 10);
}

#[test]
fn as_rng_delegates_seek() {
    use urng_core::Rng;
    let mut a = AsRng::new(Sobol64::<1>::new());
    for _ in 0..77 {
        let _ = a.nextu();
    }
    let mut b = AsRng::new(Sobol64::<1>::new());
    b.seek(77);
    assert_eq!(a, b);
    assert_eq!(a.nextu(), b.nextu());
}
