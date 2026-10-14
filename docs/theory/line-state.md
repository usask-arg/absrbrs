# The line state

The line-state layer takes one spectral line and one atmospheric condition to the four
quantities the Voigt line shape needs, q = (A, ν₀′, σ_D, y), so that the line's cross section is

σ(ν) = A · K(x, y), with x = (ν − ν₀′)/σ_D,

and K the Voigt kernel of {doc}`faddeeva`. It also gives the closed-form derivatives of q with
respect to the state and to the line parameters. The layer is crate-private until the
cross-section engine of phase 1c uses it; this page records its formulas as implemented
(`absrbrs/src/state/`).

## Inputs and units

| Input | Symbol | Unit |
|---|---|---|
| Temperature | T | K |
| Total pressure | p | Pa |
| Volume mixing ratio of the absorber | VMR | dimensionless, 0 to 1 |
| Line centre | ν₀ | cm⁻¹ |
| Intensity at T_ref | S_ref | cm⁻¹/(molecule cm⁻²) |
| Lower-state energy | E″ | cm⁻¹ |
| Air- and self-broadened half widths | γ_air, γ_self | cm⁻¹ atm⁻¹ |
| Temperature exponents of the widths | n_air, n_self | dimensionless |
| Air and self pressure shifts at T_ref | δ_air, δ_self | cm⁻¹ atm⁻¹ |
| Temperature coefficient of the air shift | δ′_air | cm⁻¹ atm⁻¹ K⁻¹ |
| Isotopologue mass | mass | g/mol |
| Partition-function ratio | ln(Q(T_ref)/Q(T)) and d ln Q/dT | dimensionless, K⁻¹ |

The line parameters are in HITRAN units at T_ref = 296 K (see {doc}`/conventions`). The
partition-function input comes from {doc}`partition-functions`, evaluated once per condition
and isotopologue. The outputs are A in m²/molecule, ν₀′ and σ_D in cm⁻¹, and y dimensionless.

## Constants

The physical constants are the exact SI values:

| Constant | Value |
|---|---|
| Speed of light c | 299 792 458 m/s |
| Boltzmann constant k_B | 1.380649 × 10⁻²³ J/K |
| Avogadro constant N_A | 6.02214076 × 10²³ mol⁻¹ |
| Planck constant h | 6.62607015 × 10⁻³⁴ J s |
| Second radiation constant c₂ = hc/k_B | 1.4387768775039338 cm K |
| Standard atmosphere | 101 325 Pa |
| Reference temperature T_ref | 296 K |

c₂ is the double nearest to hc/k_B computed from the exact constants. sasktran2 uses 1.4387769,
which is larger by 1.6 × 10⁻⁸ relative; through the Boltzmann factor that changes S(T) by
1.6 × 10⁻⁸ · c₂E″|1/T − 1/T_ref| relative, for example 1.1 × 10⁻⁷ for E″ = 3000 cm⁻¹ at
200 K.

## Formulas

Pressures enter in atm: p_atm = p/101 325, p_self = VMR · p_atm and p_air = (1 − VMR) · p_atm.

**Intensity factor.** f = S(T)/S_ref is evaluated as a single exponential times the ratio of
the stimulated-emission factors:

f = exp(ln(Q(T_ref)/Q(T)) − c₂E″(1/T − 1/T_ref)) · expm1(−c₂ν₀/T) / expm1(−c₂ν₀/T_ref).

Folding the partition-function ratio into the exponent costs one exp instead of two, and
`expm1` keeps the stimulated-emission factor accurate when c₂ν₀/T is small (microwave
lines), where 1 − exp(−c₂ν₀/T) would cancel.

**Doppler width** (the 1/e half width):

σ_D = (ν₀/c) √(2k_B T/m), with m = mass · 10⁻³/N_A kg.

**Lorentz width:**

γ_L = (T_ref/T)^n_air · γ_air · p_air + (T_ref/T)^n_self · γ_self · p_self.

**Shifted centre:**

ν₀′ = ν₀ + (δ_air + δ′_air (T − T_ref)) · p_air + δ_self · p_self.

**Shape quantities:**

A = 10⁻⁴ S_ref f/(√π σ_D) and y = γ_L/σ_D.

The factor 10⁻⁴ converts cm² to m², and 1/(√π σ_D) normalises the Voigt profile.

### Preconditions

T > 0, ν₀ > 0, mass > 0 and 0 ≤ VMR ≤ 1. The layer does not check them; its caller, the
phase 1c engine, will.

Small ν₀: `LineList` accepts any ν₀ > 0 (the readers skip records with ν₀ ≤ 0), but A and y
grow like 1/ν₀ as σ_D → 0, so q stays finite only for ν₀ ≳ 10⁻²⁸⁰ cm⁻¹; below that their true
values overflow. The ν₀ derivatives grow like 1/ν₀² and overflow likewise below about
10⁻¹⁵⁰ cm⁻¹.

## Optional fields

Three line fields are optional: n_self, δ_self and δ′_air. A NaN in an optional column means
"absent for this line". An absent n_self defaults to n_air, an absent δ_self to δ_air, and an
absent δ′_air to 0.

**The tie rule for derivatives** (refining spec §5.3):

- the derivative with respect to an absent field itself is the partial derivative taken at its
  default value;
- the derivative with respect to the field it defaults to *includes* its effect through the
  default.

So with n_self absent, ∂/∂n_air contains both the air and the self terms, and with δ_self
absent, ∂/∂δ_air contains both shifts. This keeps the Jacobian consistent with finite
differences of `with_column("n_air", …)`, which moves the defaulted values too.

## Jacobians

There are two Jacobian functions, and each also returns q:

- `line_shape_state_jacobian`: ∂q/∂T (per K), ∂q/∂p (per Pa) and ∂q/∂VMR;
- `line_shape_line_jacobian`: ∂q with respect to every line field, in `LineField` order.

They are separate so that state-only work, the common case in retrievals, never computes the
line-parameter terms. For every θ, ∂y = (∂γ_L − y ∂σ_D)/σ_D. Below, u = c₂ν₀/T,
u_ref = c₂ν₀/T_ref, g_air = (T_ref/T)^n_air γ_air and g_self = (T_ref/T)^n_self γ_self.

### State

| θ | ∂A | ∂ν₀′ | ∂σ_D | ∂γ_L |
|---|---|---|---|---|
| T | A · (−d ln Q/dT + c₂E″/T² − (u/T)/expm1(u) − 1/(2T)) | δ′_air p_air | σ_D/(2T) | −(n_air/T) g_air p_air − (n_self/T) g_self p_self |
| p | 0 | ((δ_air + δ′_air (T − T_ref))(1 − VMR) + δ_self VMR)/101 325 | 0 | (g_air (1 − VMR) + g_self VMR)/101 325 |
| VMR | 0 | p_atm (δ_self − δ_air − δ′_air (T − T_ref)) | 0 | p_atm (g_self − g_air) |

The stimulated-emission term is written as (u/T)/expm1(u), never as
exp(−u)/(1 − exp(−u)): for ultraviolet lines at low temperature expm1(u) overflows to ∞ and
the term correctly becomes 0.

### Line parameters

| θ | ∂A | ∂ν₀′ | ∂σ_D | ∂γ_L |
|---|---|---|---|---|
| ν₀ | A · ((c₂/T)/expm1(u) − (c₂/T_ref)/expm1(u_ref) − 1/ν₀) | 1 | σ_D/ν₀ | 0 |
| S_ref | 10⁻⁴ f/(√π σ_D) | 0 | 0 | 0 |
| E″ | −A c₂(1/T − 1/T_ref) | 0 | 0 | 0 |
| γ_air | 0 | 0 | 0 | (T_ref/T)^n_air p_air |
| γ_self | 0 | 0 | 0 | (T_ref/T)^n_self p_self |
| n_air | 0 | 0 | 0 | ln(T_ref/T) g_air p_air, plus ln(T_ref/T) g_self p_self if n_self is absent |
| n_self | 0 | 0 | 0 | ln(T_ref/T) g_self p_self |
| δ_air | 0 | p_air, plus p_self if δ_self is absent | 0 | 0 |
| δ′_air | 0 | (T − T_ref) p_air | 0 | 0 |
| δ_self | 0 | p_self | 0 | 0 |

∂A/∂S_ref is computed from f directly rather than as A/S_ref, so it is right for S_ref = 0.

## Validation

- **Dual numbers.** The value function is generic over its scalar type, and the tests evaluate
  it with forward-mode dual numbers carrying all 13 derivatives (T, p, VMR and the ten line
  fields). On 360 cases (8 lines, including a microwave line, a UV line whose expm1 overflows,
  E″ = n_air = S_ref = 0, and each tie on its own, at 5 temperatures from 1 K to 1500 K,
  3 pressures and 3 VMRs) every closed-form derivative agrees with the dual one to within
  10⁻¹³ (|closed| + |dual|), with the worst at 0.057 of that bound, and the values of all three
  functions are bitwise identical.
- **Arb.** `reference/data/line_state.jsonl` holds q at 192 conditions and lines, computed with
  Arb at 256 bits. The tolerance is 10⁻¹³ relative; the largest error measured is
  4.0 × 10⁻¹⁵.
