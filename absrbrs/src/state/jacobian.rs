//! Closed-form derivatives of q = (A, ν₀′, σ_D, y) with respect to the state (T, p, VMR) and to
//! every line field (spec §5.3).
//!
//! For every θ, ∂y = (∂γ_L − y ∂σ_D)/σ_D. A derivative with respect to an absent optional field
//! is the partial taken at its default; the derivative with respect to the field it defaults to
//! includes its effect through the default (∂/∂n_air includes the self term when n_self is
//! absent, and ∂/∂δ_air the self shift when δ_self is absent).

use super::{ATM, C2, Condition, FRAC_1E4_SQRT_PI, LineParameters, LineShape, T_REF, line_shape};
use crate::lines::LineField;

/// ln(Q(T_ref)/Q(T)) and d ln Q/dT at the condition's temperature.
#[allow(dead_code)] // used by the engine (phase 1c)
#[derive(Clone, Copy, Debug)]
pub(crate) struct PartitionRatio {
    pub log_q_ratio: f64,
    pub dlnq_dt: f64,
}

/// ∂q/∂T (per K), ∂q/∂p (per Pa) and ∂q/∂VMR.
#[allow(dead_code)] // used by the engine (phase 1c)
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
#[allow(dead_code)] // used by the engine (phase 1c)
pub(crate) fn line_shape_state_jacobian(
    c: &Condition<f64>,
    l: &LineParameters<f64>,
    p: PartitionRatio,
    mass: f64,
) -> (LineShape<f64>, StateJacobian) {
    let m = Intermediates::new(c, l, p, mass);
    let (q, t, vmr) = (m.q, m.t, c.vmr);
    // (u/T)/expm1(u) → 0 when expm1(u) overflows to ∞; never exp(−u)/(1 − exp(−u)).
    let dln_a_dt = -p.dlnq_dt + C2 * l.e_lower / (t * t) - m.u / t / m.u.exp_m1() - 1.0 / (2.0 * t);
    let temperature = m.derivative(
        q.amplitude * dln_a_dt,
        l.deltap_air * m.p_air,
        q.doppler_width / (2.0 * t),
        -(l.n_air / t) * m.g_air * m.p_air - (m.n_self / t) * m.g_self * m.p_self,
    );
    let pressure = m.derivative(
        0.0,
        (m.shift_air * (1.0 - vmr) + m.delta_self * vmr) / ATM,
        0.0,
        (m.g_air * (1.0 - vmr) + m.g_self * vmr) / ATM,
    );
    let vmr = m.derivative(
        0.0,
        m.p_atm * (m.delta_self - m.shift_air),
        0.0,
        m.p_atm * (m.g_self - m.g_air),
    );
    (
        q,
        StateJacobian {
            temperature,
            pressure,
            vmr,
        },
    )
}

/// q and its derivatives with respect to every line field, indexed by `LineField as usize`.
#[allow(dead_code)] // used by the engine (phase 1c)
pub(crate) fn line_shape_line_jacobian(
    c: &Condition<f64>,
    l: &LineParameters<f64>,
    p: PartitionRatio,
    mass: f64,
) -> (LineShape<f64>, [LineShape<f64>; LineField::COUNT]) {
    let m = Intermediates::new(c, l, p, mass);
    let (q, t) = (m.q, m.t);
    // The amplitude per unit S_ref, computed directly so that S_ref = 0 works. Same operations
    // as in `line_shape`.
    let stim = (-(C2 * l.nu / t)).exp_m1() / (-(C2 * l.nu / T_REF)).exp_m1();
    let boltz = (-(C2 * l.e_lower) * (1.0 / t - 1.0 / T_REF)).exp();
    let intensity = p.log_q_ratio.exp() * boltz * stim;
    let unit_amplitude = FRAC_1E4_SQRT_PI * intensity / q.doppler_width;
    let n_self_tied = l.n_self.is_none();
    let delta_self_tied = l.delta_self.is_none();

    let mut jac = [LineShape::default(); LineField::COUNT];
    for field in LineField::ALL {
        // (∂A, ∂ν₀′, ∂σ_D, ∂γ_L)
        let (da, dc, ds, dg) = match field {
            LineField::Nu => (
                // (c₂/T)/expm1(u) → 0 when expm1(u) overflows to ∞.
                q.amplitude * (C2 / t / m.u.exp_m1() - C2 / T_REF / m.u_r.exp_m1() - 1.0 / l.nu),
                1.0,
                q.doppler_width / l.nu,
                0.0,
            ),
            LineField::SRef => (unit_amplitude, 0.0, 0.0, 0.0),
            LineField::ELower => (-q.amplitude * C2 * (1.0 / t - 1.0 / T_REF), 0.0, 0.0, 0.0),
            LineField::GammaAir => (0.0, 0.0, 0.0, m.ratio.powf(l.n_air) * m.p_air),
            LineField::GammaSelf => (0.0, 0.0, 0.0, m.ratio.powf(m.n_self) * m.p_self),
            LineField::NAir => {
                let own = m.ln_ratio * m.g_air * m.p_air;
                let through_default = if n_self_tied {
                    m.ln_ratio * m.g_self * m.p_self
                } else {
                    0.0
                };
                (0.0, 0.0, 0.0, own + through_default)
            }
            LineField::NSelf => (0.0, 0.0, 0.0, m.ln_ratio * m.g_self * m.p_self),
            LineField::DeltaAir => {
                let through_default = if delta_self_tied { m.p_self } else { 0.0 };
                (0.0, m.p_air + through_default, 0.0, 0.0)
            }
            LineField::DeltapAir => (0.0, (t - T_REF) * m.p_air, 0.0, 0.0),
            LineField::DeltaSelf => (0.0, m.p_self, 0.0, 0.0),
        };
        jac[field as usize] = m.derivative(da, dc, ds, dg);
    }
    (q, jac)
}
