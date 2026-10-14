mod common;

use urng_qmc::{Sobol32, Sobol64, VanDerCorput32, VanDerCorput64};

use common::{assert_1d_stratified, assert_2d_t0, top_bits};

#[test]
fn van_der_corput_1d() {
    for (name, values) in [
        ("vdc32", top_bits(VanDerCorput32::new(), 16)),
        ("vdc64", top_bits(VanDerCorput64::new(), 16)),
    ] {
        assert_1d_stratified(&values[0], 16, name);
    }
}

#[test]
fn sobol_1d_first_64_dims() {
    for (bits, dims) in [
        (32, top_bits(Sobol32::<64>::new(), 16)),
        (64, top_bits(Sobol64::<64>::new(), 16)),
    ] {
        for (d, values) in dims.iter().enumerate() {
            assert_1d_stratified(values, 16, &format!("Sobol{bits} dimension {}", d + 1));
        }
    }
}

#[test]
fn sobol_1d_all_1024_dims() {
    for (d, values) in top_bits(Sobol64::<1024>::new(), 12).iter().enumerate() {
        assert_1d_stratified(values, 12, &format!("dimension {}", d + 1));
    }
}

#[test]
fn sobol_2d_t0() {
    let d32 = top_bits(Sobol32::<2>::new(), 14);
    assert_2d_t0(&d32[0], &d32[1], 12, "Sobol32");
    let d64 = top_bits(Sobol64::<2>::new(), 14);
    assert_2d_t0(&d64[0], &d64[1], 12, "Sobol64");
}
