use urng_qmc_core::{QuasiRng, Seek};

use crate::rd::{alpha32, alpha64};

macro_rules! weyl {
    ($name:ident, $w:ty, $golden:expr) => {
        #[doc = concat!("One-dimensional Weyl sequence `x_n = s0 + n * alpha (mod 2^", stringify!($w), "::BITS)`.")]
        ///
        /// [`golden`](Self::golden) (also the [`Default`]) uses the inverse golden ratio,
        /// the best-distributed choice. The period is `2^BITS` for any odd `alpha`, but low
        /// discrepancy needs an `alpha / 2^BITS` that is badly approximable by fractions: with
        /// `alpha = 1` the sequence is just `n`.
        ///
        /// Starts at the origin unless a seed is given (see the crate documentation).
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $name {
            state: $w,
            alpha: $w,
            seed: $w,
            index: u64,
        }

        impl $name {
            const INDEX_MASK: u64 = crate::index_mask(<$w>::BITS);

            /// Creates the sequence with step `alpha`.
            ///
            /// An even `alpha` is made odd (`alpha | 1`), since an even step shortens the
            /// period.
            #[must_use]
            pub const fn new(alpha: $w) -> Self {
                Self {
                    state: 0,
                    alpha: alpha | 1,
                    seed: 0,
                    index: 0,
                }
            }

            /// Creates the sequence whose step is the inverse golden ratio.
            #[must_use]
            pub const fn golden() -> Self {
                Self::new($golden)
            }

            /// Shifts the whole sequence by `seed`, keeping the current index.
            #[must_use]
            pub const fn with_seed(self, seed: $w) -> Self {
                Self {
                    state: seed.wrapping_add(self.alpha.wrapping_mul(self.index as $w)),
                    seed,
                    ..self
                }
            }

            /// The step.
            #[must_use]
            pub const fn alpha(&self) -> $w {
                self.alpha
            }

            /// The starting point `s0`.
            #[must_use]
            pub const fn seed(&self) -> $w {
                self.seed
            }
        }

        impl Default for $name {
            #[inline]
            fn default() -> Self {
                Self::golden()
            }
        }

        impl QuasiRng<1> for $name {
            type Word = $w;

            #[inline]
            fn next_point(&mut self) -> [$w; 1] {
                let x = self.state;
                self.state = self.state.wrapping_add(self.alpha);
                self.index = self.index.wrapping_add(1) & Self::INDEX_MASK;
                [x]
            }
        }

        impl Seek for $name {
            #[inline]
            fn seek(&mut self, index: u64) {
                let index = index & Self::INDEX_MASK;
                self.state = self.seed.wrapping_add(self.alpha.wrapping_mul(index as $w));
                self.index = index;
            }

            #[inline]
            fn index(&self) -> u64 {
                self.index
            }
        }
    };
}

weyl!(Weyl32, u32, alpha32::<1>()[0]);
weyl!(Weyl64, u64, alpha64::<1>()[0]);
