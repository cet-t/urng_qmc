use urng_qmc_core::{QuasiRng, Seek};

/// Largest dimension [`Halton32`] and [`Halton64`] support.
pub const HALTON_MAX_DIM: usize = 32;

/// The base of each coordinate: the first 32 primes.
const PRIMES: [u64; HALTON_MAX_DIM] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131,
];

/// Most base-`b` digits of a `u64` for `b >= 3`: `3^40 < 2^64 < 3^41`.
const MAX_DIGITS: usize = 41;

/// `SCALE[j][k] = floor((2^128 - 1) / PRIMES[j]^(k + 1))` for the odd bases. Row 0 (base 2)
/// is unused: that coordinate is a bit reversal.
///
/// A `static`, not a `const`, so that runtime indexing reads it in place.
static SCALE: [[u128; MAX_DIGITS]; HALTON_MAX_DIM] = scale_table();

const fn scale_table() -> [[u128; MAX_DIGITS]; HALTON_MAX_DIM] {
    let mut t = [[0; MAX_DIGITS]; HALTON_MAX_DIM];
    let mut j = 1;
    while j < HALTON_MAX_DIM {
        let mut s = u128::MAX;
        let mut k = 0;
        while k < MAX_DIGITS {
            s /= PRIMES[j] as u128;
            t[j][k] = s;
            k += 1;
        }
        j += 1;
    }
    t
}

/// The radical inverse of `n` in base `PRIMES[j]` (`j >= 1`) as a 128-bit fixed-point
/// number.
///
/// Digit `d_k` adds `d_k * floor((2^128 - 1) / b^k)`, so the sum never overflows and falls
/// short of the exact value by at most `(b - 1)` per digit: about `2^-118` in all. The top
/// 32 or 64 bits are therefore the exact value rounded down, except within that distance
/// below a multiple of the output resolution.
///
/// The base is a const parameter so that the digit loop divides by a constant (a multiply),
/// which makes it several times faster than dividing by a runtime base.
#[inline]
fn radical_inverse_in<const J: usize>(mut n: u64) -> u128 {
    let b = PRIMES[J];
    let scale = &SCALE[J];
    let mut acc = 0;
    let mut k = 0;
    while n != 0 {
        acc += (n % b) as u128 * scale[k];
        n /= b;
        k += 1;
    }
    acc
}

/// [`radical_inverse_in`] for coordinate `j` in `1..HALTON_MAX_DIM`.
#[inline]
fn radical_inverse(j: usize, n: u64) -> u128 {
    macro_rules! dispatch {
        ($($j:literal)*) => {
            match j {
                $($j => radical_inverse_in::<$j>(n),)*
                _ => unreachable!(),
            }
        };
    }
    dispatch!(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31)
}

macro_rules! halton {
    ($name:ident, $w:ty) => {
        #[doc = concat!("Halton sequence with `", stringify!($w), "` coordinates.")]
        ///
        /// Coordinate `j` (0-based) is the radical inverse of the index in the `j + 1`-th
        /// prime (2, 3, 5, ..., 131). Coordinate 0 is exactly the van der Corput sequence; the
        /// others are within one unit of the exact value rounded down.
        ///
        /// Coordinates with large bases are strongly correlated over the first points, so
        /// prefer Sobol beyond about 8 dimensions. `D` must be in `1..=HALTON_MAX_DIM`;
        /// anything else fails to compile:
        ///
        /// ```compile_fail
        #[doc = concat!("let _ = urng_qmc::", stringify!($name), "::<33>::new();")]
        /// ```
        ///
        /// The first point is the origin (see the crate documentation). The index wraps
        /// modulo `2^BITS` like the other sequences, although Halton itself has no period.
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $name<const D: usize> {
            index: u64,
        }

        impl<const D: usize> $name<D> {
            const INDEX_MASK: u64 = crate::index_mask(<$w>::BITS);

            /// Creates the sequence at index 0.
            #[must_use]
            pub const fn new() -> Self {
                const {
                    assert!(
                        D >= 1 && D <= HALTON_MAX_DIM,
                        "Halton supports 1 to 32 dimensions"
                    )
                };
                Self { index: 0 }
            }

            #[inline]
            fn point(n: u64) -> [$w; D] {
                core::array::from_fn(|j| {
                    if j == 0 {
                        (n as $w).reverse_bits()
                    } else {
                        (radical_inverse(j, n) >> (128 - <$w>::BITS)) as $w
                    }
                })
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
                let p = Self::point(self.index);
                self.index = self.index.wrapping_add(1) & Self::INDEX_MASK;
                p
            }
        }

        impl<const D: usize> Seek for $name<D> {
            #[inline]
            fn seek(&mut self, index: u64) {
                self.index = index & Self::INDEX_MASK;
            }

            #[inline]
            fn index(&self) -> u64 {
                self.index
            }
        }
    };
}

halton!(Halton32, u32);
halton!(Halton64, u64);
