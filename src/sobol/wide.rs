//! Sobol with all 21201 Joe–Kuo dimensions (`sobol-21201` feature).
//!
//! [`Sobol32`](super::Sobol32) builds one direction table per `D` at compile time. Past a
//! couple of thousand dimensions that single evaluation trips the deny-by-default
//! `long_running_const_eval` lint, and the lint is checked in the crate that names the
//! type, so users would have to allow it. Here the table is split into chunks of 256
//! dimensions, each its own constant and well under the limit, and a point is updated chunk
//! by chunk.

use urng_qmc_core::{DigitalNet, QuasiRng, Seek};

use super::table_21201::{MAX_DIM, PARAMS};
use super::{directions32, directions64};

/// Largest dimension [`SobolWide32`] and [`SobolWide64`] support.
pub const SOBOL_WIDE_MAX_DIM: usize = MAX_DIM;

/// Dimensions per chunk.
const CHUNK: usize = 256;

macro_rules! chunk_builder {
    ($name:ident, $w:ty, $bits:literal, $directions:ident) => {
        /// Direction numbers `v[k][j]` of dimensions `c * CHUNK + j`; zero past `MAX_DIM`.
        const fn $name(c: usize) -> [[$w; CHUNK]; $bits] {
            let mut v = [[0; CHUNK]; $bits];
            let mut j = 0;
            while j < CHUNK {
                let d = c * CHUNK + j;
                if d < MAX_DIM {
                    let column = $directions(if d == 0 { None } else { Some(PARAMS[d - 1]) });
                    let mut k = 0;
                    while k < $bits {
                        v[k][j] = column[k];
                        k += 1;
                    }
                }
                j += 1;
            }
            v
        }
    };
}

chunk_builder!(chunk32, u32, 32, directions32);
chunk_builder!(chunk64, u64, 64, directions64);

/// One constant per chunk, so that each is evaluated separately (well under the lint's
/// limit) and only when a type with enough dimensions is used.
struct Chunk<const C: usize>;

impl<const C: usize> Chunk<C> {
    const V32: &'static [[u32; CHUNK]; 32] = &chunk32(C);
    const V64: &'static [[u64; CHUNK]; 64] = &chunk64(C);
}

// Behind a call, a chunk is only evaluated when the call runs: constants named directly in
// an initializer are evaluated even on branches that are not taken.
const fn chunk32_ref<const C: usize>() -> &'static [[u32; CHUNK]; 32] {
    Chunk::<C>::V32
}

const fn chunk64_ref<const C: usize>() -> &'static [[u64; CHUNK]; 64] {
    Chunk::<C>::V64
}

/// Stands in for the chunks past the last dimension a type uses.
const UNUSED32: &[[u32; CHUNK]; 32] = &[[0; CHUNK]; 32];
const UNUSED64: &[[u64; CHUNK]; 64] = &[[0; CHUNK]; 64];

/// The first `$used` chunks, then `$unused`.
macro_rules! chunks {
    ($get:ident, $unused:ident, $used:expr; $($c:literal)*) => {
        [$(if $c < $used { $get::<$c>() } else { $unused }),*]
    };
}

const N_CHUNKS: usize = MAX_DIM.div_ceil(CHUNK);

macro_rules! sobol_wide {
    ($name:ident, $w:ty, $bits:literal, $get:ident, $unused:ident) => {
        #[doc = concat!("Sobol sequence with `", stringify!($w), "` coordinates and up to all 21201 Joe–Kuo dimensions.")]
        ///
        /// The same points as the narrower Sobol type for `D <= 1024`, which is lighter there.
        /// The direction numbers are built at compile time in chunks of 256 dimensions, only
        /// for the chunks `D` needs, and kept in read-only data (`BITS` words per dimension:
        /// 10.4 MiB for `SobolWide64<21201>`). Building all of them takes on the order of a
        /// minute per width. Prefer [`fill_point`](QuasiRng::fill_point) to returning large
        /// points by value.
        ///
        /// `D` must be in `1..=SOBOL_WIDE_MAX_DIM`; anything else fails to compile:
        ///
        /// ```compile_fail
        #[doc = concat!("let _ = urng_qmc::", stringify!($name), "::<21202>::new();")]
        /// ```
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $name<const D: usize> {
            x: [$w; D],
            index: u64,
        }

        impl<const D: usize> $name<D> {
            /// The chunks covering dimensions `0..D`. Evaluated per `D`, when the type is
            /// used, rather than with the crate.
            const CHUNKS: [&'static [[$w; CHUNK]; $bits]; N_CHUNKS] =
                chunks!($get, $unused, D.div_ceil(CHUNK); 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82);

            /// Creates the sequence at index 0.
            #[must_use]
            pub const fn new() -> Self {
                const {
                    assert!(
                        D >= 1 && D <= MAX_DIM,
                        "SobolWide supports 1 to 21201 dimensions"
                    )
                };
                Self { x: [0; D], index: 0 }
            }

            /// XORs direction number `k` of every dimension into `x`.
            #[inline]
            fn xor_row(x: &mut [$w; D], k: usize) {
                for (xs, chunk) in x.chunks_mut(CHUNK).zip(&Self::CHUNKS) {
                    for (x, v) in xs.iter_mut().zip(&chunk[k]) {
                        *x ^= *v;
                    }
                }
            }

            #[inline]
            fn advance(&mut self) {
                let c = self.index.trailing_ones() as usize;
                if c >= $bits {
                    self.x = [0; D];
                    self.index = 0;
                } else {
                    Self::xor_row(&mut self.x, c);
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
                    Self::xor_row(&mut self.x, gray.trailing_zeros() as usize);
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

sobol_wide!(SobolWide32, u32, 32, chunk32_ref, UNUSED32);
sobol_wide!(SobolWide64, u64, 64, chunk64_ref, UNUSED64);
