//! The line-state layer (spec §4.2, §5): from a condition and a line to the shape quantities
//! q = (A, ν₀′, σ_D, y), with closed-form derivatives in `jacobian`.
//!
//! Pressures are in atm: p_atm = p/101325, p_self = VMR·p_atm, p_air = (1 − VMR)·p_atm.
//!
//! - Intensity factor: f = S(T)/S_ref = exp(ln(Q(T_ref)/Q(T)) − c₂E″(1/T − 1/T_ref)) ·
//!   expm1(−c₂ν₀/T)/expm1(−c₂ν₀/T_ref).
//! - Doppler width: σ_D = (ν₀/c) √(2k_B T/m), with m = mass·10⁻³/N_A kg.
//! - Lorentz width: γ_L = (T_ref/T)^{n_air} γ_air p_air + (T_ref/T)^{n_self} γ_self p_self.
//! - Shifted centre: ν₀′ = ν₀ + (δ_air + δ′_air (T − T_ref)) p_air + δ_self p_self.
//! - Shape quantities: A = 10⁻⁴ S_ref f/(√π σ_D) and y = γ_L/σ_D, so that the cross section is
//!   σ(ν) = A·K(x, y) in m²/molecule.
//!
//! An absent n_self defaults to n_air, and an absent δ_self to δ_air.

// Nothing outside the tests calls into this module yet; test builds still lint dead code.
#![cfg_attr(not(test), allow(dead_code, unused_imports))] // used by the engine (phase 1c)

#[cfg(test)]
mod dual;
mod jacobian;
mod real;

pub(crate) use jacobian::{
    PartitionRatio, StateJacobian, line_shape_line_jacobian, line_shape_state_jacobian,
};
pub(crate) use real::Real;

use crate::lines::{LineField, LineList};

/// Reference temperature of the line parameters, K.
pub(crate) const T_REF: f64 = 296.0;
/// One standard atmosphere, Pa.
pub(crate) const ATM: f64 = 101_325.0;
/// Second radiation constant hc/k_B, cm K (exact SI constants; nearest double).
pub(crate) const C2: f64 = 1.438_776_877_503_933_8;
/// Speed of light, m/s.
pub(crate) const C_LIGHT: f64 = 299_792_458.0;
/// Boltzmann constant, J/K.
pub(crate) const K_B: f64 = 1.380_649e-23;
/// Avogadro constant, mol⁻¹.
pub(crate) const N_A: f64 = 6.022_140_76e23;
/// 10⁻⁴/√π: the cm² to m² conversion and the Voigt normalisation.
const FRAC_1E4_SQRT_PI: f64 = 1e-4 * crate::faddeeva::FRAC_1_SQRT_PI;

/// An atmospheric condition.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Condition<R> {
    /// Temperature T, K.
    pub temperature: R,
    /// Total pressure p, Pa.
    pub pressure: R,
    /// Volume mixing ratio of the absorber (dimensionless, 0 to 1).
    pub vmr: R,
}

/// The parameters of one line, in HITRAN units at T_ref = 296 K.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LineParameters<R> {
    /// Line centre ν₀, cm⁻¹.
    pub nu: R,
    /// Intensity S_ref at T_ref, cm⁻¹/(molecule cm⁻²).
    pub s_ref: R,
    /// Lower-state energy E″, cm⁻¹.
    pub e_lower: R,
    /// Air-broadened half width γ_air, cm⁻¹ atm⁻¹.
    pub gamma_air: R,
    /// Self-broadened half width γ_self, cm⁻¹ atm⁻¹.
    pub gamma_self: R,
    /// Temperature exponent n_air of γ_air (dimensionless).
    pub n_air: R,
    /// Temperature exponent n_self of γ_self (dimensionless); `None` when absent, and then
    /// n_air is used.
    pub n_self: Option<R>,
    /// Air pressure shift δ_air at T_ref, cm⁻¹ atm⁻¹.
    pub delta_air: R,
    /// Temperature coefficient δ′_air of the air shift, cm⁻¹ atm⁻¹ K⁻¹ (0 when absent).
    pub deltap_air: R,
    /// Self pressure shift δ_self, cm⁻¹ atm⁻¹; `None` when absent, and then δ_air is used.
    pub delta_self: Option<R>,
}

impl LineParameters<f64> {
    /// The parameters of line `i` of a list: n_self and δ_self are `None` when absent, and an
    /// absent δ′_air is 0.
    ///
    /// # Panics
    ///
    /// If `i >= list.len()`.
    pub(crate) fn from_list(list: &LineList, i: usize) -> LineParameters<f64> {
        let optional = |field| list.has_value(field, i).then(|| list.value(field, i));
        LineParameters {
            nu: list.value(LineField::Nu, i),
            s_ref: list.value(LineField::SRef, i),
            e_lower: list.value(LineField::ELower, i),
            gamma_air: list.value(LineField::GammaAir, i),
            gamma_self: list.value(LineField::GammaSelf, i),
            n_air: list.value(LineField::NAir, i),
            n_self: optional(LineField::NSelf),
            delta_air: list.value(LineField::DeltaAir, i),
            deltap_air: list.value(LineField::DeltapAir, i),
            delta_self: optional(LineField::DeltaSelf),
        }
    }
}

/// The shape quantities q of one line at one condition: σ(ν) = A·K(x, y) with
/// x = (ν − ν₀′)/σ_D. Also used for derivatives of q, component by component.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct LineShape<R> {
    /// A, m²/molecule.
    pub amplitude: R,
    /// The pressure-shifted centre ν₀′, cm⁻¹.
    pub centre: R,
    /// The Doppler 1/e half width σ_D, cm⁻¹.
    pub doppler_width: R,
    /// y = γ_L/σ_D.
    pub y: R,
}

/// The intensity factor f = S(T)/S_ref (see the module docs). Shared by [`line_shape`] and the
/// S_ref derivative, so the two agree bitwise.
fn intensity_factor<R: Real>(t: R, nu: R, e_lower: R, log_q_ratio: R) -> R {
    let f = R::from_f64;
    let tr = f(T_REF);
    let stim = (-(f(C2) * nu / t)).exp_m1() / (-(f(C2) * nu / tr)).exp_m1();
    (log_q_ratio - f(C2) * e_lower * (f(1.0) / t - f(1.0) / tr)).exp() * stim
}

/// q for one line at one condition (see the module docs for the formulas). `log_q_ratio` is
/// ln(Q(T_ref)/Q(T)) and `mass` the isotopologue mass in g/mol.
///
/// Preconditions: T > 0, ν₀ > 0, mass > 0 and 0 ≤ VMR ≤ 1 (the engine checks them). q stays
/// finite for ν₀ ≳ 1e-280 cm⁻¹, although [`LineList`] accepts any ν₀ > 0: A and y grow like
/// 1/ν₀ (σ_D → 0), so below that their true values overflow. (The ν₀ derivatives, which grow
/// like 1/ν₀², overflow likewise below about 1e-150 cm⁻¹.)
pub(crate) fn line_shape<R: Real>(
    c: &Condition<R>,
    l: &LineParameters<R>,
    log_q_ratio: R,
    mass: f64,
) -> LineShape<R> {
    let f = R::from_f64;
    let p = c.pressure / f(ATM);
    let p_self = c.vmr * p;
    let p_air = (f(1.0) - c.vmr) * p;
    let t = c.temperature;
    let tr = f(T_REF);
    let intensity = intensity_factor(t, l.nu, l.e_lower, log_q_ratio);
    // 2k_B/m with m = mass·10⁻³/N_A kg.
    let sigma_d = l.nu / f(C_LIGHT) * (f(2.0 * K_B * N_A * 1e3 / mass) * t).sqrt();
    let ratio = tr / t;
    let n_self = l.n_self.unwrap_or(l.n_air);
    let delta_self = l.delta_self.unwrap_or(l.delta_air);
    let gamma =
        ratio.powf(l.n_air) * l.gamma_air * p_air + ratio.powf(n_self) * l.gamma_self * p_self;
    LineShape {
        amplitude: f(FRAC_1E4_SQRT_PI) * l.s_ref * intensity / sigma_d,
        centre: l.nu + (l.delta_air + l.deltap_air * (t - tr)) * p_air + delta_self * p_self,
        doppler_width: sigma_d,
        y: gamma / sigma_d,
    }
}

#[cfg(test)]
mod dual_checks {
    use super::dual::Dual;
    use super::*;

    /// T, p, VMR, then the line fields in `LineField::ALL` order.
    const N: usize = 3 + LineField::COUNT;
    type D = Dual<N>;

    /// (ν₀, S_ref, E″, γ_air, γ_self, n_air, n_self, δ_air, δ′_air, δ_self).
    #[rustfmt::skip]
    type Line = (f64, f64, f64, f64, f64, f64, Option<f64>, f64, f64, Option<f64>);

    #[rustfmt::skip]
    const LINES: [Line; 8] = [
        (13100.5, 7.4e-24, 1200.0, 0.045, 0.047, 0.72, Some(0.70), -0.008, 1e-5, Some(-0.007)),
        (6300.25, 1.6e-23, 300.0, 0.07, 0.09, 0.75, None, -0.004, 0.0, None),
        // A microwave line: c₂ν₀/T is small.
        (0.75, 1e-24, 45.0, 0.1, 0.4, 0.68, None, 0.0, 0.0, None),
        (2140.0, 4e-19, 8000.0, 0.05, 0.06, 0.7, Some(0.6), -0.002, 0.0, Some(0.003)),
        // Only δ_self tied, then only n_self tied: the tie rule field by field.
        (4300.0, 2e-21, 500.0, 0.06, 0.08, 0.73, Some(0.66), -0.005, 2e-5, None),
        (1600.0, 3e-20, 150.0, 0.08, 0.35, 0.69, None, -0.003, 0.0, Some(-0.006)),
        // E″ = 0, n_air = 0 and S_ref = 0 (A = 0, but ∂A/∂S_ref is not).
        (950.0, 0.0, 0.0, 0.06, 0.07, 0.0, None, -0.001, 0.0, None),
        // A UV line: expm1(c₂ν₀/T) overflows to ∞ at low T.
        (40000.0, 5e-22, 100.0, 0.05, 0.055, 0.75, Some(0.7), -0.01, 3e-6, Some(-0.012)),
    ];
    /// Isotopologue masses, g/mol, one per line.
    const MASSES: [f64; 8] = [32.0, 44.0, 18.0, 28.0, 64.0, 16.0, 46.0, 48.0];
    const TEMPERATURES: [f64; 5] = [1.0, 50.0, 150.0, 296.0, 1500.0];
    const PRESSURES: [f64; 3] = [0.0, 1000.0, 101_325.0];
    const VMRS: [f64; 3] = [0.0, 0.01, 1.0];
    const LOG_Q_RATIO: f64 = 0.37;
    const DLNQ_DT: f64 = 0.0042;

    fn var(v: f64, field: LineField) -> D {
        D::variable(v, 3 + field as usize)
    }

    fn f64_line(line: &Line) -> LineParameters<f64> {
        let &(nu, s_ref, e_lower, gamma_air, gamma_self, n_air, n_self, delta_air, dp, ds) = line;
        LineParameters {
            nu,
            s_ref,
            e_lower,
            gamma_air,
            gamma_self,
            n_air,
            n_self,
            delta_air,
            deltap_air: dp,
            delta_self: ds,
        }
    }

    /// Dual parameters; with `seed_defaults`, absent n_self and δ_self become present at their
    /// defaults (n_air, δ_air), seeded as their own variables.
    fn dual_line(l: &LineParameters<f64>, seed_defaults: bool) -> LineParameters<D> {
        let optional = |v: Option<f64>, default: f64, field| match v {
            Some(v) => Some(var(v, field)),
            None if seed_defaults => Some(var(default, field)),
            None => None,
        };
        LineParameters {
            nu: var(l.nu, LineField::Nu),
            s_ref: var(l.s_ref, LineField::SRef),
            e_lower: var(l.e_lower, LineField::ELower),
            gamma_air: var(l.gamma_air, LineField::GammaAir),
            gamma_self: var(l.gamma_self, LineField::GammaSelf),
            n_air: var(l.n_air, LineField::NAir),
            n_self: optional(l.n_self, l.n_air, LineField::NSelf),
            delta_air: var(l.delta_air, LineField::DeltaAir),
            deltap_air: var(l.deltap_air, LineField::DeltapAir),
            delta_self: optional(l.delta_self, l.delta_air, LineField::DeltaSelf),
        }
    }

    fn components<R: Copy>(q: &LineShape<R>) -> [(&'static str, R); 4] {
        [
            ("A", q.amplitude),
            ("centre", q.centre),
            ("sigma_D", q.doppler_width),
            ("y", q.y),
        ]
    }

    /// Relative difference, scaled so that the bound is |closed − dual| ≤ 1e-13·(|closed| +
    /// |dual|) + 1e-300.
    fn check(closed: f64, dual: f64, what: &str, worst: &mut f64) {
        let err = (closed - dual).abs() / ((closed.abs() + dual.abs()) * 1e-13 + 1e-300);
        *worst = worst.max(if err.is_nan() { f64::INFINITY } else { err });
        assert!(
            (closed - dual).abs() <= 1e-13 * (closed.abs() + dual.abs()) + 1e-300,
            "{what}: closed {closed:e}, dual {dual:e}"
        );
    }

    #[test]
    fn closed_forms_match_dual_numbers() {
        let mut worst = 0.0f64; // in units of the 1e-13 relative bound
        let mut cases = 0;
        for (line, mass) in LINES.iter().zip(MASSES) {
            let l = f64_line(line);
            let tied = l.n_self.is_none() || l.delta_self.is_none();
            for t in TEMPERATURES {
                for p in PRESSURES {
                    for vmr in VMRS {
                        cases += 1;
                        let c = Condition {
                            temperature: t,
                            pressure: p,
                            vmr,
                        };
                        let cd = Condition {
                            temperature: D::variable(t, 0),
                            pressure: D::variable(p, 1),
                            vmr: D::variable(vmr, 2),
                        };
                        let mut lq = D::constant(LOG_Q_RATIO);
                        lq.d[0] = -DLNQ_DT;
                        let pr = PartitionRatio {
                            log_q_ratio: LOG_Q_RATIO,
                            dlnq_dt: DLNQ_DT,
                        };
                        let q = line_shape(&c, &l, LOG_Q_RATIO, mass);
                        let qd = line_shape(&cd, &dual_line(&l, false), lq, mass);
                        let (qs, js): (_, StateJacobian) =
                            line_shape_state_jacobian(&c, &l, pr, mass);
                        let (ql, jl) = line_shape_line_jacobian(&c, &l, pr, mass);
                        let name = format!("ν₀ = {}, T = {t}, p = {p}, VMR = {vmr}", l.nu);
                        for ((what, v), ((_, vd), ((_, vs), (_, vl)))) in
                            components(&q).into_iter().zip(
                                components(&qd)
                                    .into_iter()
                                    .zip(components(&qs).into_iter().zip(components(&ql))),
                            )
                        {
                            assert_eq!(v.to_bits(), vd.v.to_bits(), "{name}: {what} dual");
                            assert_eq!(v.to_bits(), vs.to_bits(), "{name}: {what} state");
                            assert_eq!(v.to_bits(), vl.to_bits(), "{name}: {what} line");
                        }
                        let state = [("T", js.temperature), ("p", js.pressure), ("VMR", js.vmr)];
                        for (k, (theta, jac)) in state.into_iter().enumerate() {
                            for ((what, cl), (_, du)) in
                                components(&jac).into_iter().zip(components(&qd))
                            {
                                check(
                                    cl,
                                    du.d[k],
                                    &format!("{name}: ∂{what}/∂{theta}"),
                                    &mut worst,
                                );
                            }
                        }
                        let qt = tied.then(|| line_shape(&cd, &dual_line(&l, true), lq, mass));
                        for field in LineField::ALL {
                            let k = 3 + field as usize;
                            // The partial at the default needs the run with the default seeded.
                            let at_default = match field {
                                LineField::NSelf => l.n_self.is_none(),
                                LineField::DeltaSelf => l.delta_self.is_none(),
                                _ => false,
                            };
                            let reference = if at_default {
                                qt.as_ref().unwrap()
                            } else {
                                &qd
                            };
                            for ((what, cl), (_, du)) in components(&jl[field as usize])
                                .into_iter()
                                .zip(components(reference))
                            {
                                let theta = field.name();
                                check(
                                    cl,
                                    du.d[k],
                                    &format!("{name}: ∂{what}/∂{theta}"),
                                    &mut worst,
                                );
                            }
                        }
                    }
                }
            }
        }
        eprintln!(
            "line state: {cases} cases; worst |closed − dual| = {worst:.3} × the 1e-13 bound"
        );
    }

    /// `from_list` reads absent optional fields as `None` (n_self, δ_self) or 0 (δ′_air).
    #[test]
    fn from_list_marks_absent_fields() {
        use crate::lines::LineColumns;
        let list = LineList::from_columns(LineColumns {
            molecule: vec![2, 2],
            isotopologue: vec![1, 1],
            nu: vec![1000.0, 1001.0],
            s_ref: vec![1e-20, 2e-20],
            e_lower: vec![10.0, 20.0],
            gamma_air: vec![0.07, 0.08],
            gamma_self: vec![0.09, 0.1],
            n_air: vec![0.7, 0.75],
            n_self: Some(vec![0.65, f64::NAN]),
            delta_air: vec![-0.001, -0.002],
            deltap_air: Some(vec![f64::NAN, 2e-5]),
            delta_self: None,
        })
        .unwrap();
        let a = LineParameters::from_list(&list, 0);
        assert_eq!(
            (a.nu, a.n_self, a.deltap_air, a.delta_self),
            (1000.0, Some(0.65), 0.0, None)
        );
        let b = LineParameters::from_list(&list, 1);
        assert_eq!(
            (b.nu, b.n_self, b.deltap_air, b.delta_self),
            (1001.0, None, 2e-5, None)
        );
        assert_eq!(
            (b.s_ref, b.e_lower, b.gamma_air, b.gamma_self),
            (2e-20, 20.0, 0.08, 0.1)
        );
        assert_eq!((b.n_air, b.delta_air), (0.75, -0.002));
    }
}

#[cfg(test)]
mod reference {
    use super::*;
    use crate::test_data::reference_file;

    /// q against Arb at 256 bits (`reference/data/line_state.jsonl`).
    #[test]
    fn matches_arb_reference() {
        let Some(text) = reference_file("line_state.jsonl") else {
            return;
        };
        let mut worst = [0.0f64; 4];
        let mut n = 0;
        for line in text.lines().filter(|l| !l.trim().is_empty()) {
            let r: serde_json::Value = serde_json::from_str(line).unwrap();
            // A missing key is an error; only `null` means "absent" (n_self, δ_self).
            let input = |name: &str| match r["inputs"].get(name) {
                Some(serde_json::Value::Null) => None,
                Some(v) => Some(v.as_f64().unwrap_or_else(|| panic!("input {name}: {v}"))),
                None => panic!("{}: no input {name}", r["case"]),
            };
            let get = |name: &str| input(name).unwrap_or_else(|| panic!("input {name} is null"));
            let c = Condition {
                temperature: get("temperature"),
                pressure: get("pressure"),
                vmr: get("vmr"),
            };
            let l = LineParameters {
                nu: get("nu"),
                s_ref: get("s_ref"),
                e_lower: get("e_lower"),
                gamma_air: get("gamma_air"),
                gamma_self: get("gamma_self"),
                n_air: get("n_air"),
                n_self: input("n_self"), // null when absent
                delta_air: get("delta_air"),
                deltap_air: get("deltap_air"),
                delta_self: input("delta_self"),
            };
            let q = line_shape(&c, &l, get("log_q_ratio"), get("mass"));
            let actual = [q.amplitude, q.centre, q.doppler_width, q.y];
            let names = ["amplitude", "centre", "doppler_width", "y"];
            for (k, (name, v)) in names.into_iter().zip(actual).enumerate() {
                let out = &r["outputs"][name];
                let hi = out["hi"][0].as_f64().unwrap();
                let lo = out["lo"][0].as_f64().unwrap();
                // An exact zero (y at p = 0) must be reproduced exactly.
                let err = if hi == 0.0 {
                    if v == 0.0 { 0.0 } else { f64::INFINITY }
                } else {
                    ((v - hi) - lo).abs() / hi.abs()
                };
                worst[k] = if err.is_nan() {
                    f64::INFINITY
                } else {
                    worst[k].max(err)
                };
                assert!(
                    err <= 1e-13,
                    "{}: {name} = {v:e}, Arb {hi:e}: {err:e}",
                    r["case"]
                );
            }
            n += 1;
        }
        assert_eq!(n, 192);
        eprintln!(
            "line state vs Arb, {n} records: worst relative error A {:e}, ν₀′ {:e}, σ_D {:e}, y {:e}",
            worst[0], worst[1], worst[2], worst[3]
        );
    }
}
