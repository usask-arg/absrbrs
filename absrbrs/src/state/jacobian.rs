//! Closed-form derivatives of q = (A, ν₀′, σ_D, y) with respect to the state (T, p, VMR) and to
//! every line field (spec §5.3).
//!
//! For every θ, ∂y = (∂γ_L − y ∂σ_D)/σ_D. A derivative with respect to an absent optional field
//! is the partial taken at its default; the derivative with respect to the field it defaults to
//! includes its effect through the default (∂/∂n_air includes the self term when n_self is
//! absent, and ∂/∂δ_air the self shift when δ_self is absent).

use super::{
    ATM, C2, Condition, FRAC_1E4_SQRT_PI, LineParameters, LineShape, T_REF, intensity_factor,
    line_shape,
};
use crate::lines::LineField;

/// The partition-function input of one isotopologue at the condition's temperature T.
#[derive(Clone, Copy, Debug)]
pub(crate) struct PartitionRatio {
    /// ln(Q(T_ref)/Q(T)) (dimensionless).
    pub log_q_ratio: f64,
    /// d ln Q/dT at T, K⁻¹, from the same interpolant as Q (so ∂/∂T of `log_q_ratio` is
    /// −`dlnq_dt`).
    pub dlnq_dt: f64,
}

/// ∂q/∂T (per K), ∂q/∂p (per Pa) and ∂q/∂VMR.
#[derive(Clone, Copy, Debug)]
pub(crate) struct StateJacobian {
    pub temperature: LineShape<f64>,
    pub pressure: LineShape<f64>,
    pub vmr: LineShape<f64>,
}

/// The intermediates shared by the two Jacobians.
struct Intermediates {
    q: LineShape<f64>,
    t: f64,
    p_atm: f64,
    p_air: f64,
    p_self: f64,
    ratio: f64,
    ln_ratio: f64,
    n_self: f64,
    delta_self: f64,
    g_air: f64,
    g_self: f64,
    shift_air: f64,
    u: f64,
    u_r: f64,
}

impl Intermediates {
    fn new(c: &Condition<f64>, l: &LineParameters<f64>, p: PartitionRatio, mass: f64) -> Self {
        let q = line_shape(c, l, p.log_q_ratio, mass);
        let t = c.temperature;
        let p_atm = c.pressure / ATM;
        let ratio = T_REF / t;
        let n_self = l.n_self.unwrap_or(l.n_air);
        Intermediates {
            q,
            t,
            p_atm,
            p_air: (1.0 - c.vmr) * p_atm,
            p_self: c.vmr * p_atm,
            ratio,
            ln_ratio: ratio.ln(),
            n_self,
            delta_self: l.delta_self.unwrap_or(l.delta_air),
            g_air: ratio.powf(l.n_air) * l.gamma_air,
            g_self: ratio.powf(n_self) * l.gamma_self,
            shift_air: l.delta_air + l.deltap_air * (t - T_REF),
            u: C2 * l.nu / t,
            u_r: C2 * l.nu / T_REF,
        }
    }

    /// ∂q from ∂A, ∂ν₀′, ∂σ_D and ∂γ_L.
    fn derivative(
        &self,
        amplitude: f64,
        centre: f64,
        doppler_width: f64,
        gamma: f64,
    ) -> LineShape<f64> {
        LineShape {
            amplitude,
            centre,
            doppler_width,
            y: (gamma - self.q.y * doppler_width) / self.q.doppler_width,
        }
    }
}

/// q and its derivatives with respect to T, p and VMR.
pub(crate) fn line_shape_state_jacobian(
    c: &Condition<f64>,
    l: &LineParameters<f64>,
    p: PartitionRatio,
    mass: f64,
) -> (LineShape<f64>, StateJacobian) {
    let im = Intermediates::new(c, l, p, mass);
    let (q, t, x) = (im.q, im.t, c.vmr);
    // (u/T)/expm1(u) → 0 when expm1(u) overflows to ∞; never exp(−u)/(1 − exp(−u)).
    let dln_a_dt =
        -p.dlnq_dt + C2 * l.e_lower / (t * t) - im.u / t / im.u.exp_m1() - 1.0 / (2.0 * t);
    let d_temperature = im.derivative(
        q.amplitude * dln_a_dt,
        l.deltap_air * im.p_air,
        q.doppler_width / (2.0 * t),
        -(l.n_air / t) * im.g_air * im.p_air - (im.n_self / t) * im.g_self * im.p_self,
    );
    let d_pressure = im.derivative(
        0.0,
        (im.shift_air * (1.0 - x) + im.delta_self * x) / ATM,
        0.0,
        (im.g_air * (1.0 - x) + im.g_self * x) / ATM,
    );
    let d_vmr = im.derivative(
        0.0,
        im.p_atm * (im.delta_self - im.shift_air),
        0.0,
        im.p_atm * (im.g_self - im.g_air),
    );
    (
        q,
        StateJacobian {
            temperature: d_temperature,
            pressure: d_pressure,
            vmr: d_vmr,
        },
    )
}

/// q and its derivatives with respect to every line field, indexed by `LineField as usize`.
pub(crate) fn line_shape_line_jacobian(
    c: &Condition<f64>,
    l: &LineParameters<f64>,
    p: PartitionRatio,
    mass: f64,
) -> (LineShape<f64>, [LineShape<f64>; LineField::COUNT]) {
    let im = Intermediates::new(c, l, p, mass);
    let (q, t) = (im.q, im.t);
    // The amplitude per unit S_ref, computed directly so that S_ref = 0 works.
    let intensity = intensity_factor(t, l.nu, l.e_lower, p.log_q_ratio);
    let unit_amplitude = FRAC_1E4_SQRT_PI * intensity / q.doppler_width;
    let n_self_tied = l.n_self.is_none();
    let delta_self_tied = l.delta_self.is_none();

    let mut jac = [LineShape::default(); LineField::COUNT];
    for field in LineField::ALL {
        // (∂A, ∂ν₀′, ∂σ_D, ∂γ_L)
        let (da, dc, ds, dg) = match field {
            LineField::Nu => (
                // (c₂/T)/expm1(u) → 0 when expm1(u) overflows to ∞.
                q.amplitude * (C2 / t / im.u.exp_m1() - C2 / T_REF / im.u_r.exp_m1() - 1.0 / l.nu),
                1.0,
                q.doppler_width / l.nu,
                0.0,
            ),
            LineField::SRef => (unit_amplitude, 0.0, 0.0, 0.0),
            LineField::ELower => (-q.amplitude * C2 * (1.0 / t - 1.0 / T_REF), 0.0, 0.0, 0.0),
            LineField::GammaAir => (0.0, 0.0, 0.0, im.ratio.powf(l.n_air) * im.p_air),
            LineField::GammaSelf => (0.0, 0.0, 0.0, im.ratio.powf(im.n_self) * im.p_self),
            LineField::NAir => {
                let own = im.ln_ratio * im.g_air * im.p_air;
                let through_default = if n_self_tied {
                    im.ln_ratio * im.g_self * im.p_self
                } else {
                    0.0
                };
                (0.0, 0.0, 0.0, own + through_default)
            }
            LineField::NSelf => (0.0, 0.0, 0.0, im.ln_ratio * im.g_self * im.p_self),
            LineField::DeltaAir => {
                let through_default = if delta_self_tied { im.p_self } else { 0.0 };
                (0.0, im.p_air + through_default, 0.0, 0.0)
            }
            LineField::DeltapAir => (0.0, (t - T_REF) * im.p_air, 0.0, 0.0),
            LineField::DeltaSelf => (0.0, im.p_self, 0.0, 0.0),
        };
        jac[field as usize] = im.derivative(da, dc, ds, dg);
    }
    (q, jac)
}
