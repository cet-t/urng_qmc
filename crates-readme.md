# urng_qmc

Quasi-random (low-discrepancy) sequences on top of [`urng_core`](https://crates.io/crates/urng_core). `no_std`, no `alloc`, no `unsafe`.

The traits come from [`urng_qmc_core`](https://crates.io/crates/urng_qmc_core) and are re-exported here, so this crate is the only dependency you need.

| Type | Dimensions | Notes |
| --- | --- | --- |
| `Weyl32`, `Weyl64` | 1 | `x_n = s0 + n * alpha`. `golden()` uses the inverse golden ratio |
| `Rd32<D>`, `Rd64<D>` | 1 to 16 | Roberts' R_d: one Weyl sequence per coordinate |
| `VanDerCorput32`, `VanDerCorput64` | 1 | Base-2 radical inverse |
| `Sobol32<D>`, `Sobol64<D>` | 1 to 1024 | Joe–Kuo direction numbers, Gray-code order |
| `SobolWide32<D>`, `SobolWide64<D>` | 1 to 21201 | Same, all Joe–Kuo dimensions; `sobol-21201` feature |
| `Halton32<D>`, `Halton64<D>` | 1 to 32 | Radical inverses in the first 32 primes |
| `Shifted<Q, D>` | any | Random shift (Cranley–Patterson), for randomized QMC |
| `DigitalShifted<Q, D>` | any | Random XOR (digital shift), for digital nets |

Every type implements `QuasiRng<D>` and `Seek`. The van der Corput and Sobol types and `DigitalShifted` also implement `DigitalNet<D>`.

## Usage

```rust
use urng_qmc::{QuasiRng, Sobol64};

// Estimate π with 2^16 points.
let mut q = Sobol64::<2>::new();
let mut inside = 0u32;
for _ in 0..(1 << 16) {
    let [x, y] = q.next_point_f(); // [f64; 2] in [0, 1)
    if x * x + y * y < 1.0 {
        inside += 1;
    }
}
let pi = 4.0 * inside as f64 / (1 << 16) as f64;
assert!((pi - core::f64::consts::PI).abs() < 1e-3);
```

`next_point` returns the raw fixed-point words; `next_point_f` and `next_point_in` convert them through `urng_core::Word::to_float` and `randf`.

Randomize a sequence with any `urng_core::Rng`, to estimate the error from independent copies:

```rust
use urng_core::Rng;
use urng_qmc::{DigitalShifted, QuasiRng, Sobol64};

struct SplitMix(u64);
impl Rng for SplitMix {
    type Word = u64;
    fn nextu(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let z = (self.0 ^ (self.0 >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        let z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

let mut rng = SplitMix(42);
let estimates: Vec<f64> = (0..8)
    .map(|_| {
        let mut q = DigitalShifted::new(Sobol64::<2>::new(), &mut rng);
        (0..1024).map(|_| { let [x, y] = q.next_point_f(); x * y }).sum::<f64>() / 1024.0
    })
    .collect();
let mean = estimates.iter().sum::<f64>() / 8.0; // ∫∫ xy = 1/4
assert!((mean - 0.25).abs() < 1e-3);
```

A one-dimensional sequence becomes a `urng_core::Rng` through `AsRng`:

```rust
use urng_core::Rng;
use urng_qmc::{AsRng, Sobol32};

let mut rng = AsRng::new(Sobol32::<1>::new());
let _: f32 = rng.nextf();
```

## Notes

- **Not suitable for cryptographic use.**
- **Draw whole points, not successive values.** Successive values of a low-discrepancy sequence are strongly correlated. Code that draws several values per sample from one `Rng` (2-D points, Box–Muller, rejection sampling) breaks silently with `AsRng`: estimating π from pairs of successive van der Corput values gives about 3.29. Use a `D`-dimensional sequence instead.
- **The first point is the origin** for Sobol, van der Corput, and unseeded Weyl and R_d. It is kept on purpose: Sobol is balanced over blocks of `2^m` points, and dropping the first point breaks that (Owen, *On dropping the first Sobol' point*, 2020). `ln(x)` or an inverse CDF can return infinity there; scramble the sequence, or call `seek(1)` if you must skip it. Even scrambled, coordinates below the float resolution convert to `0.0`.
- **Use the 64-bit types for many points.** `to_float` keeps 23 bits for `f32`, so more than `2^23` points of a 32-bit sequence start to collide after conversion.
- Every sequence has period `2^BITS`, and `seek` reduces its index modulo `2^BITS`.
- Sobol's direction table is built at compile time for the `D` you use and placed in read-only data: `D * BITS` words (512 KiB for `Sobol64<1024>`). `SobolWide` builds it in chunks of 256 dimensions: about 8 s of compile time for `SobolWide64<3000>`, about a minute per width for all 21201 dimensions. Enabling the feature without using the type costs nothing.
- Halton coordinates with large bases are strongly correlated over the first points; prefer Sobol beyond about 8 dimensions.
