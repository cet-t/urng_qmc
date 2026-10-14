use urng_qmc_core::{DigitalNet, QuasiRng, Seek};

macro_rules! vdc {
    ($name:ident, $w:ty) => {
        #[doc = concat!("Base-2 van der Corput sequence with `", stringify!($w), "` values: `x_n = (n as ", stringify!($w), ").reverse_bits()`.")]
        ///
        /// Holds the same point set as the first coordinate of Sobol, in a different order.
        /// The first point is the origin (see the crate documentation).
        #[derive(Clone, Debug, PartialEq, Eq, Default)]
        pub struct $name {
            index: u64,
        }

        impl $name {
            const INDEX_MASK: u64 = crate::index_mask(<$w>::BITS);

            /// Creates the sequence at index 0.
            #[must_use]
            pub const fn new() -> Self {
                Self { index: 0 }
            }
        }

        impl QuasiRng<1> for $name {
            type Word = $w;

            #[inline]
            fn next_point(&mut self) -> [$w; 1] {
                let x = (self.index as $w).reverse_bits();
                self.index = self.index.wrapping_add(1) & Self::INDEX_MASK;
                [x]
            }
        }

        impl Seek for $name {
            #[inline]
            fn seek(&mut self, index: u64) {
                self.index = index & Self::INDEX_MASK;
            }

            #[inline]
            fn index(&self) -> u64 {
                self.index
            }
        }

        impl DigitalNet<1> for $name {}
    };
}

vdc!(VanDerCorput32, u32);
vdc!(VanDerCorput64, u64);
