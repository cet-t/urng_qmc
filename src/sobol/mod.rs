mod table;
#[cfg(feature = "sobol-21201")]
mod table_21201;
#[cfg(feature = "sobol-21201")]
mod wide;

use urng_qmc_core::{DigitalNet, QuasiRng, Seek};

use self::table::{MAX_DIM, PARAMS};
#[cfg(feature = "sobol-21201")]
pub use self::wide::{SOBOL_WIDE_MAX_DIM, SobolWide32, SobolWide64};

/// Largest dimension [`Sobol32`] and [`Sobol64`] support. The `sobol-21201` feature adds
/// `SobolWide32` and `SobolWide64` for up to 21201.
pub const SOBOL_MAX_DIM: usize = MAX_DIM;

macro_rules! directions {
    ($name:ident, $w:ty, $bits:literal) => {
        /// The direction numbers `v_0..v_{BITS-1}` of one dimension: `row` holds its
        /// Joe–Kuo parameters `(s, a, m)`, or is `None` for dimension 1.
        const fn $name<const N: usize>(row: Option<(u8, u32, [u32; N])>) -> [$w; $bits] {
            let mut v = [0; $bits];
            let mut k = 0;
            match row {
                None => {
                    while k < $bits {
                        v[k] = 1 << ($bits - 1 - k);
                        k += 1;
                    }
                }
                Some((s, a, m)) => {
                    let s = s as usize;
                    while k < $bits {
                        v[k] = if k < s {
                            (m[k] as $w) << ($bits - 1 - k)
                        } else {
                            let mut x = v[k - s] ^ (v[k - s] >> s);
                            // Bit `s - 1 - i` of `a` is the coefficient a_i, i in 1..s.
                            // Visiting only the set bits halves the const-evaluation work.
                            let mut bits = a;
                            while bits != 0 {
                                x ^= v[k - (s - 1 - bits.trailing_zeros() as usize)];
                                bits &= bits - 1;
                            }
                            x
                        };
                        k += 1;
                    }
                }
            }
            v
        }
    };
}

directions!(directions32, u32, 32);
directions!(directions64, u64, 64);

macro_rules! sobol {
    ($name:ident, $w:ty, $bits:literal, $build:ident, $directions:ident) => {
        /// Direction numbers `v[k][d]`: bit position `k` outermost, so that one update
        /// reads a contiguous row.
        const fn $build<const D: usize>() -> [[$w; D]; $bits] {
            assert!(D >= 1 && D <= MAX_DIM, "Sobol supports 1 to 1024 dimensions");
            let mut v = [[0; D]; $bits];
            let mut d = 0;
            while d < D {
                let column = $directions(if d == 0 { None } else { Some(PARAMS[d - 1]) });
                let mut k = 0;
                while k < $bits {
                    v[k][d] = column[k];
                    k += 1;
                }
                d += 1;
            }
            v
        }

        #[doc = concat!("Sobol sequence with `", stringify!($w), "` coordinates and the Joe–Kuo `new-joe-kuo-6.21201` direction numbers.")]
        ///
        /// Points come in Gray-code order, like SciPy's: the first `2^m` points are the same
        /// set as in natural order, so every block of `2^m` points starting at a multiple of
        /// `2^m` is a `(t, m, D)`-net. One update costs `O(D)` and a [`Seek`] `O(D * BITS)`.
        ///
        /// `D` must be in `1..=SOBOL_MAX_DIM` (1024); anything else fails to compile. The
        /// `sobol-21201` feature adds `SobolWide32` / `SobolWide64` for more dimensions.
        ///
        /// ```compile_fail
        #[doc = concat!("let _ = urng_qmc::", stringify!($name), "::<0>::new();")]
        /// ```
        ///
        /// The first point is the origin (see the crate documentation).
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $name<const D: usize> {
            x: [$w; D],
            index: u64,
        }

        impl<const D: usize> $name<D> {
            // A reference, not a value: a `const` array is copied wherever it is used,
            // which in debug builds means the whole table on every index.
            const V: &'static [[$w; D]; $bits] = &$build::<D>();

            /// Creates the sequence at index 0.
            #[must_use]
            pub const fn new() -> Self {
                const { assert!(D >= 1 && D <= MAX_DIM, "Sobol supports 1 to 1024 dimensions") };
                Self { x: [0; D], index: 0 }
            }

            #[inline]
            fn advance(&mut self) {
                let c = self.index.trailing_ones() as usize;
                if c >= $bits {
                    // All 2^BITS points used: start over.
                    self.x = [0; D];
                    self.index = 0;
                } else {
                    for (x, v) in self.x.iter_mut().zip(&Self::V[c]) {
                        *x ^= *v;
                    }
                    self.index += 1;
                }
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
                let out = self.x;
                self.advance();
                out
            }

            #[inline]
            fn fill_point(&mut self, out: &mut [$w; D]) {
                *out = self.x;
                self.advance();
            }
        }

        impl<const D: usize> Seek for $name<D> {
            fn seek(&mut self, index: u64) {
                let index = index & crate::index_mask($bits);
                self.x = [0; D];
                let mut gray = index ^ (index >> 1);
                while gray != 0 {
                    let k = gray.trailing_zeros() as usize;
                    for (x, v) in self.x.iter_mut().zip(&Self::V[k]) {
                        *x ^= *v;
                    }
                    gray &= gray - 1;
                }
                self.index = index;
            }

            #[inline]
            fn index(&self) -> u64 {
                self.index
            }
        }

        impl<const D: usize> DigitalNet<D> for $name<D> {}
    };
}

sobol!(Sobol32, u32, 32, build32, directions32);
sobol!(Sobol64, u64, 64, build64, directions64);
