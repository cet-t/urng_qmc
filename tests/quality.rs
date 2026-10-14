//! "Works as a quasi-random sequence": discrepancy and integration error.

mod common;

use std::f64::consts::PI;

use urng_qmc::{DigitalShifted, Halton64, QuasiRng, Rd64, Shifted, Sobol64};

use common::{SplitMix64, l2_star, l2_star_random, slope, units};

fn discrepancy_ratio<Q: QuasiRng<D>, const D: usize>(mut q: Q, n: usize) -> f64 {
    let points: Vec<[f64; D]> = (0..n).map(|_| units(q.next_point())).collect();
    l2_star(&points) / l2_star_random(D, n)
}

/// Thresholds from section 8 of docs/spec.md: measured value plus about 25 %.
#[test]
fn l2_star_discrepancy() {
    let n = 1024;
    let cases = [
        ("Sobol D=2", discrepancy_ratio(Sobol64::<2>::new(), n), 0.25),
        ("R_d D=2", discrepancy_ratio(Rd64::<2>::new(), n), 0.25),
        ("Sobol D=4", discrepancy_ratio(Sobol64::<4>::new(), n), 0.25),
        ("R_d D=4", discrepancy_ratio(Rd64::<4>::new(), n), 0.95),
        ("Sobol D=8", discrepancy_ratio(Sobol64::<8>::new(), n), 0.85),
        ("R_d D=8", discrepancy_ratio(Rd64::<8>::new(), n), 1.03),
        (
            "Halton D=2",
            discrepancy_ratio(Halton64::<2>::new(), n),
            0.25,
        ),
        (
            "Halton D=4",
            discrepancy_ratio(Halton64::<4>::new(), n),
            0.35,
        ),
        (
            "Halton D=8",
            discrepancy_ratio(Halton64::<8>::new(), n),
            1.37,
        ),
    ];
    for (name, ratio, limit) in cases {
        eprintln!("{name}: {ratio:.3} of random (limit {limit})");
        assert!(ratio <= limit, "{name}: {ratio:.3} > {limit}");
    }
}

/// `prod (pi / 2) sin(pi x_j)`, whose integral over the unit cube is 1.
fn integrand<const D: usize>(p: [f64; D]) -> f64 {
    p.iter().map(|x| PI / 2.0 * (PI * x).sin()).product()
}

/// Slope of log2(RMSE) against log2(N) for N = 2^8..=2^16, over 16 randomized copies
/// made by `scramble`.
fn convergence_slope<Q, const D: usize>(scramble: impl Fn(&mut SplitMix64) -> Q) -> f64
where
    Q: QuasiRng<D>,
{
    let ns: Vec<usize> = (8..=16).map(|m| 1 << m).collect();
    let mut sq_err = vec![0.0; ns.len()];
    let mut rng = SplitMix64(2026);
    let reps = 16;
    for _ in 0..reps {
        let mut q = scramble(&mut rng);
        let mut sum = 0.0;
        let mut next = 0;
        for i in 1..=*ns.last().unwrap() {
            sum += integrand(units(q.next_point()));
            if i == ns[next] {
                sq_err[next] += (sum / i as f64 - 1.0).powi(2);
                next += 1;
            }
        }
    }
    let xs: Vec<f64> = ns.iter().map(|&n| (n as f64).log2()).collect();
    let ys: Vec<f64> = sq_err
        .iter()
        .map(|e| (e / reps as f64).sqrt().log2())
        .collect();
    slope(&xs, &ys)
}

#[test]
fn integration_error_converges_faster_than_monte_carlo() {
    let sobol = convergence_slope(|rng| DigitalShifted::new(Sobol64::<2>::new(), rng));
    let rd = convergence_slope(|rng| Shifted::new(Rd64::<2>::new(), rng));
    eprintln!("RMSE slope: Sobol {sobol:.3}, R_d {rd:.3} (Monte Carlo is about -0.5)");
    assert!(sobol < -0.8, "Sobol slope {sobol:.3}");
    assert!(rd < -0.8, "R_d slope {rd:.3}");
}
