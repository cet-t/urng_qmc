use urng_core::Rng;
use urng_qmc_core::{DigitalNet, QmcWord, QuasiRng, Seek};

/// Cranley–Patterson rotation: adds the same random word to coordinate `j` of every
/// point, modulo `2^BITS`.
///
/// Every point becomes uniform on `[0, 1)^D` while the spacing between points is kept, so
/// averages over independent rotations give an unbiased estimate and its error
/// (randomized quasi-Monte Carlo). It suits Weyl and R_d, whose lattice structure it keeps.
/// For digital nets such as Sobol, use [`DigitalShifted`] instead: a rotation breaks the
/// `(t, m, s)`-net property.
///
/// ```
/// use urng_qmc::{QuasiRng, Rd64, Shifted};
/// # struct SplitMix(u64);
/// # impl urng_core::Rng for SplitMix {
/// #     type Word = u64;
/// #     fn nextu(&mut self) -> u64 {
/// #         self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
/// #         let z = (self.0 ^ (self.0 >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
/// #         let z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
/// #         z ^ (z >> 31)
/// #     }
/// # }
///
/// let mut rng = SplitMix(1);
/// let mut q = Shifted::new(Rd64::<2>::new(), &mut rng);
/// let [x, y] = q.next_point_f();
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Shifted<Q, const D: usize>
where
    Q: QuasiRng<D>,
{
    inner: Q,
    shift: [Q::Word; D],
}

/// Digital shift: XORs the same random word into coordinate `j` of every point.
///
/// Like [`Shifted`] this makes every point uniform, but it keeps the `(t, m, s)`-net
/// property, so it accepts only [`DigitalNet`]s and is itself a [`DigitalNet`].
///
/// ```
/// use urng_qmc::{DigitalShifted, QuasiRng, Sobol64};
/// # struct Counter(u32);
/// # impl urng_core::Rng for Counter {
/// #     type Word = u32;
/// #     fn nextu(&mut self) -> u32 {
/// #         self.0 = self.0.wrapping_mul(0x9E37_79B9).wrapping_add(1);
/// #         self.0
/// #     }
/// # }
///
/// // A 32-bit generator can scramble a 64-bit sequence: two draws per word.
/// let mut rng = Counter(7);
/// let mut q = DigitalShifted::new(Sobol64::<3>::new(), &mut rng);
/// let _: [f64; 3] = q.next_point_f();
/// ```
///
/// Sequences that are not digital nets are rejected at compile time:
///
/// ```compile_fail
/// use urng_qmc::{DigitalShifted, Rd64};
/// let _ = DigitalShifted::with_mask(Rd64::<2>::new(), [1, 2]);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DigitalShifted<Q, const D: usize>
where
    Q: DigitalNet<D>,
{
    inner: Q,
    mask: [Q::Word; D],
}

macro_rules! wrapper_common {
    ($name:ident, $bound:ident, $field:ident, $op:expr, $what:literal) => {
        impl<Q: $bound<D>, const D: usize> $name<Q, D> {
            #[doc = concat!("Wraps `inner` with a ", $what, " drawn from `rng`, one word per coordinate.")]
            ///
            /// The generator may have either word width (see
            /// [`QmcWord::from_rng`](urng_qmc_core::QmcWord::from_rng)).
            #[must_use]
            pub fn new<R: Rng + ?Sized>(inner: Q, rng: &mut R) -> Self
            where
                R::Word: QmcWord,
            {
                let $field = core::array::from_fn(|_| Q::Word::from_rng(rng));
                Self { inner, $field }
            }

            #[doc = concat!("The ", $what, " applied to coordinate `j`.")]
            #[must_use]
            pub const fn $field(&self) -> &[Q::Word; D] {
                &self.$field
            }

            /// Returns a reference to the wrapped sequence.
            #[must_use]
            pub const fn get_ref(&self) -> &Q {
                &self.inner
            }

            /// Returns the wrapped sequence.
            #[must_use]
            pub fn into_inner(self) -> Q {
                self.inner
            }
        }

        impl<Q: $bound<D>, const D: usize> QuasiRng<D> for $name<Q, D> {
            type Word = Q::Word;

            #[inline]
            fn next_point(&mut self) -> [Q::Word; D] {
                let mut p = self.inner.next_point();
                for (x, &s) in p.iter_mut().zip(&self.$field) {
                    *x = $op(*x, s);
                }
                p
            }

            #[inline]
            fn fill_point(&mut self, out: &mut [Q::Word; D]) {
                self.inner.fill_point(out);
                for (x, &s) in out.iter_mut().zip(&self.$field) {
                    *x = $op(*x, s);
                }
            }
        }

        /// Delegates to the wrapped sequence: the transform is the same for every point.
        impl<Q: $bound<D> + Seek, const D: usize> Seek for $name<Q, D> {
            #[inline]
            fn seek(&mut self, index: u64) {
                self.inner.seek(index);
            }

            #[inline]
            fn index(&self) -> u64 {
                self.inner.index()
            }
        }
    };
}

wrapper_common!(Shifted, QuasiRng, shift, QmcWord::wrapping_add, "shift");
wrapper_common!(DigitalShifted, DigitalNet, mask, |x, m| x ^ m, "XOR mask");

impl<Q: QuasiRng<D>, const D: usize> Shifted<Q, D> {
    /// Wraps `inner` with a given shift.
    #[must_use]
    pub const fn with_shift(inner: Q, shift: [Q::Word; D]) -> Self {
        Self { inner, shift }
    }
}

impl<Q: DigitalNet<D>, const D: usize> DigitalShifted<Q, D> {
    /// Wraps `inner` with a given XOR mask.
    #[must_use]
    pub const fn with_mask(inner: Q, mask: [Q::Word; D]) -> Self {
        Self { inner, mask }
    }
}

impl<Q: DigitalNet<D>, const D: usize> DigitalNet<D> for DigitalShifted<Q, D> {}
