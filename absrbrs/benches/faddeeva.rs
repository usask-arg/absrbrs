//! Throughput of the Faddeeva kernels: one point in each Fast-tier region, and the grid points
//! of representative line windows (the phase 0 workload: O₂ A-band and CO₂ 1.6 µm lines at the
//! surface, near 30 km and near 65 km, 0.005 cm⁻¹ grid, 25 cm⁻¹ window), against sasktran2's
//! 8-pole J-pole approximation.
//!
//! The 8-pole J-pole baseline below is sasktran2's, copied here for comparison only. It is
//! licensed as follows:
//!
//! MIT License
//!
//! Copyright (c) 2024 University of Saskatchewan Atmospheric Research Group
//!
//! Permission is hereby granted, free of charge, to any person obtaining a copy of
//! this software and associated documentation files (the “Software”), to deal in
//! the Software without restriction, including without limitation the rights to
//! use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
//! of the Software, and to permit persons to whom the Software is furnished to do
//! so, subject to the following conditions:
//!
//! The above copyright notice and this permission notice shall be included in all
//! copies or substantial portions of the Software.
//!
//! THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//! SOFTWARE.

use std::hint::black_box;

use absrbrs::faddeeva::{Complex64, Tier, kernel, kernel_into};
use criterion::{Criterion, Throughput, criterion_group, criterion_main};

/// sasktran2's J-pole constants (rust/sasktran2-rs/src/math/errorfunctions/optimized.rs):
/// Copyright (c) 2024 University of Saskatchewan Atmospheric Research Group, MIT licence.
const JPOLE_B: [Complex64; 4] = [
    Complex64::new(0.00383968430671409, -0.0119854387180615),
    Complex64::new(-0.321597857664957, -0.218883985607935),
    Complex64::new(2.55515264319988, 0.613958600684469),
    Complex64::new(-2.73739446984183, 5.69007914897806),
];
const JPOLE_C: [Complex64; 4] = [
    Complex64::new(2.51506776338386, -1.60713668042405),
    Complex64::new(-1.68985621846204, -1.66471695485661),
    Complex64::new(0.981465428659098, -1.70017951305004),
    Complex64::new(-0.322078795578047, -1.71891780447016),
];
const FRAC_1_SQRT_PI: f64 = 0.564_189_583_547_756_3;

/// The J-pole kernel triplet: the same 8-pole sum as sasktran2's, evaluated in plain complex
/// arithmetic. The other four poles are −conj of these, with conj coefficients. sasktran2's own
/// derivative path precomputes per-line pole-pair coefficients and uses f64x8 SIMD under its
/// nightly-only feature, which this baseline does not represent (an informal scratch measurement,
/// not reproducible from the repository, put a scalar pole-paired version at 12.9 ns per point
/// against this one's 13.4 ns).
fn jpole_kernel(x: f64, y: f64) -> (f64, f64, f64) {
    let z = Complex64::new(x, y);
    let (mut s, mut ds) = (Complex64::new(0.0, 0.0), Complex64::new(0.0, 0.0));
    for (b, c) in JPOLE_B.iter().zip(JPOLE_C.iter()) {
        for (b, c) in [(*b, *c), (b.conj(), -c.conj())] {
            let inv = (z - c).inv();
            s += b * inv;
            ds += b * inv * inv;
        }
    }
    let i = Complex64::new(0.0, 1.0);
    let w = -i * FRAC_1_SQRT_PI * s;
    let dw = i * FRAC_1_SQRT_PI * ds;
    (w.re, dw.re, -dw.im)
}

const K_B: f64 = 1.380649e-23;
const N_A: f64 = 6.02214076e23;
const C: f64 = 299_792_458.0;

/// (name, centre cm⁻¹, mass g/mol, γ_air cm⁻¹/atm, n_air).
const LINES: [(&str, f64, f64, f64, f64); 2] = [
    ("o2", 13100.0, 31.98983, 0.045, 0.72),
    ("co2", 6300.0, 43.98983, 0.070, 0.75),
];
/// (name, temperature K, pressure Pa).
const CONDITIONS: [(&str, f64, f64); 3] = [
    ("surface", 250.0, 101325.0),
    ("30km", 225.0, 1000.0),
    ("65km", 200.0, 10.0),
];

/// x at every grid point of a line's window, and the line's y.
fn window(centre: f64, mass: f64, gamma: f64, n: f64, t: f64, p: f64) -> (Vec<f64>, f64) {
    let sigma_d = centre / C * (2.0 * K_B * t * N_A * 1e3 / mass).sqrt();
    let y = gamma * (296.0 / t).powf(n) * p / 101325.0 / sigma_d;
    (
        (-5000..=5000)
            .map(|i| f64::from(i) * 0.005 / sigma_d)
            .collect(),
        y,
    )
}

fn bench_regions(c: &mut Criterion) {
    let mut group = c.benchmark_group("point");
    let points = [
        ("weideman", 2.0, 0.5),
        ("taylor", 3.0, 1e-3),
        ("exp_cf5", 10.0, 1e-3),
        ("cf4", 8.5, 0.5),
        ("cf3", 12.0, 1.0),
        ("cf2", 30.0, 1.0),
        ("cf1", 300.0, 1.0),
        ("cf0", 3000.0, 1.0),
    ];
    for (name, x, y) in points {
        group.bench_function(format!("fast/{name}"), |b| {
            b.iter(|| kernel(black_box(x), black_box(y), Tier::Fast))
        });
        group.bench_function(format!("accurate/{name}"), |b| {
            b.iter(|| kernel(black_box(x), black_box(y), Tier::Accurate))
        });
        group.bench_function(format!("jpole/{name}"), |b| {
            b.iter(|| jpole_kernel(black_box(x), black_box(y)))
        });
    }
    group.finish();
}

fn bench_lines(c: &mut Criterion) {
    let mut group = c.benchmark_group("line");
    for (line, centre, mass, gamma, n) in LINES {
        for (cond, t, p) in CONDITIONS {
            let (x, y) = window(centre, mass, gamma, n, t, p);
            let len = x.len();
            let (mut k, mut kx, mut ky) = (vec![0.0; len], vec![0.0; len], vec![0.0; len]);
            group.throughput(Throughput::Elements(len as u64));
            for (tier_name, tier) in [("fast", Tier::Fast), ("accurate", Tier::Accurate)] {
                group.bench_function(format!("{tier_name}/{line}_{cond}"), |b| {
                    b.iter(|| {
                        kernel_into(black_box(&x), y, tier, &mut k, &mut kx, &mut ky).unwrap();
                        black_box((&mut k, &mut kx, &mut ky));
                    })
                });
            }
            group.bench_function(format!("jpole/{line}_{cond}"), |b| {
                b.iter(|| {
                    for (i, &xi) in x.iter().enumerate() {
                        (k[i], kx[i], ky[i]) = jpole_kernel(black_box(xi), y);
                    }
                    black_box((&mut k, &mut kx, &mut ky));
                })
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_regions, bench_lines);
criterion_main!(benches);
