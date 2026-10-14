#![doc = include_str!("../crates-readme.md")]
#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod rd;
mod sobol;
mod vdc;
mod weyl;

pub use urng_qmc_core::*;

pub use crate::rd::{RD_MAX_DIM, Rd32, Rd64};
pub use crate::sobol::{SOBOL_MAX_DIM, Sobol32, Sobol64};
pub use crate::vdc::{VanDerCorput32, VanDerCorput64};
pub use crate::weyl::{Weyl32, Weyl64};

/// Mask that keeps an index in `[0, 2^bits)`.
const fn index_mask(bits: u32) -> u64 {
    u64::MAX >> (64 - bits)
}
