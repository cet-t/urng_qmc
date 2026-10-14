mod table;

use urng_qmc_core::{QuasiRng, Seek};

use self::table::{ALPHA, MAX_DIM};

/// Largest dimension [`Rd32`] and [`Rd64`] support.
pub const RD_MAX_DIM: usize = MAX_DIM;

/// The 64-bit steps for dimension `D`.
pub(crate) const fn alpha64<const D: usize>() -> [u64; D] {
    assert!(D >= 1 && D <= MAX_DIM, "R_d supports 1 to 16 dimensions");
    let mut a = [0; D];
    let mut j = 0;
    while j < D {
        a[j] = ALPHA[D - 1][j];
        j += 1;
    }
    a
}

/// The 32-bit steps for dimension `D`: the high halves of the 64-bit ones, made odd.
pub(crate) const fn alpha32<const D: usize>() -> [u32; D] {
    let wide = alpha64::<D>();
    let mut a = [0; D];
    let mut j = 0;
    while j < D {
        a[j] = (wide[j] >> 32) as u32 | 1;
        j += 1;
    }
    a
}

macro_rules! rd {
    ($name:ident, $w:ty, $alpha:ident) => {
        #[doc = concat!("Roberts' R_d sequence with `", stringify!($w), "` coordinates.")]
        ///
        /// Coordinate `j` (1-based) is a Weyl sequence with step `floor(phi_D^-j * 2^BITS) | 1`,
        /// where `phi_D` is the positive root of `x^(D+1) = x + 1`. `D` must be in
        /// `1..=RD_MAX_DIM`; anything else fails to compile:
        ///
        /// ```compile_fail
        #[doc = concat!("let _ = urng_qmc::", stringify!($name), "::<17>::new();")]
        /// ```
        ///
        /// Starts at the origin unless a seed is given (see the crate documentation).
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $name<const D: usize> {
            state: [$w; D],
            seed: [$w; D],
            index: u64,
        }

        impl<const D: usize> $name<D> {
            /// The step of each coordinate.
            pub const ALPHA: [$w; D] = $alpha::<D>();

            const INDEX_MASK: u64 = crate::index_mask(<$w>::BITS);

            /// Creates the sequence, starting at the origin.
            #[must_use]
            pub const fn new() -> Self {
                const { assert!(D >= 1 && D <= MAX_DIM, "R_d supports 1 to 16 dimensions") };
                Self {
                    state: [0; D],
                    seed: [0; D],
                    index: 0,
                }
            }

            /// Shifts coordinate `j` of every point by `seed[j]`, keeping the current index.
            #[must_use]
            pub fn with_seed(mut self, seed: [$w; D]) -> Self {
                self.seed = seed;
                self.seek(self.index);
                self
            }

            /// The starting point.
            #[must_use]
            pub const fn seed(&self) -> &[$w; D] {
                &self.seed
            }
        }

        impl<const D: usize> Default for $name<D> {
            #[inline]
            fn default() -> Self {
                Self::new()
            }
        }

        impl<const D: usize> QuasiRng<D> for $name<D> {
            type Word = $w;

            #[inline]
            fn next_point(&mut self) -> [$w; D] {
                let out = self.state;
                for (x, a) in self.state.iter_mut().zip(Self::ALPHA) {
                    *x = x.wrapping_add(a);
                }
                self.index = self.index.wrapping_add(1) & Self::INDEX_MASK;
                out
            }
        }

        impl<const D: usize> Seek for $name<D> {
            #[inline]
            fn seek(&mut self, index: u64) {
                let index = index & Self::INDEX_MASK;
                for ((x, s), a) in self.state.iter_mut().zip(self.seed).zip(Self::ALPHA) {
                    *x = s.wrapping_add(a.wrapping_mul(index as $w));
                }
                self.index = index;
            }

            #[inline]
            fn index(&self) -> u64 {
                self.index
            }
        }
    };
}

rd!(Rd32, u32, alpha32);
rd!(Rd64, u64, alpha64);
