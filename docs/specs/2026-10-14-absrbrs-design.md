# absrbrs — design specification

- **Date:** 2026-10-14
- **Status:** approved design, pre-implementation
- **Scope of this spec:** the overall architecture and roadmap, plus the full design of phase 1 (line-by-line Voigt cross sections with Jacobians)

**Later decisions.** Decisions taken after 2026-10-14 are recorded in `docs/plans/phase-0-followups.md`, which supersedes this spec where they differ: §5.5 (TIPS table grid and the comparison with HAPI), §6.1 (the kernel returns K, K_x and K_y; in the Taylor region these are the formula's partial derivatives), §6.3 (the Fast-tier algorithm is chosen) and §8.6 (HITRAN snippets are blocked; loader tests use synthetic records).

## 1. Purpose

absrbrs is a Rust library with a Python interface for molecular absorption cross sections.
Its name reads as both "absorbers" and "absrb-rs".

Its distinguishing features are:

1. **First-class linearization.** Every cross section has a Jacobian with respect to the atmospheric state (temperature, pressure, volume mixing ratio) and, on request, with respect to the spectroscopic line parameters.
   - State derivatives are the primary use case. They feed the sasktran2 weighting-function framework and must be cheap.
   - Line-parameter derivatives are opt-in. They make it possible to fit lab spectra, and they must never slow down the state-only path.
2. **Speed** on atmospheric radiative-transfer workloads.
3. **Bounded, validated accuracy.** The line-shape kernels are validated against arbitrary-precision references (Arb), with accuracy maps enforced in CI.

absrbrs will consolidate the absorption code now in sasktran2 (`rust/sasktran2-rs/src/optical/`), and sasktran2 will eventually depend on it as a Rust crate. Over time the scope grows towards a radis-like package: line mixing, non-Voigt profiles, fast algorithms for very large line lists, continua and collision-induced absorption.

### 1.1 Phase 1 includes

- An SoA line list with HITRAN `.par` and HAPI `.data`/`.header` readers.
- Built-in TIPS-2021 partition functions, with an analytic d ln Q/dT and no Python callbacks.
- HITRAN-convention temperature and pressure scaling, with air and self broadening and optional n_self, δ_self and δ′_air columns.
- A Faddeeva function with derivatives at two accuracy tiers, Fast and Accurate.
- A windowed Voigt broadening engine with configurable cutoffs and optional pedestal subtraction.
- Forward-mode Jacobians for state parameters and user-selected line parameters.
- Adjoint (vector–Jacobian) gradients for all state and line parameters.
- A Python package with a functional API, a reusable calculator object and a `LineList` class.
- Validation infrastructure, benchmarks against sasktran2 and HAPI, CI, docs and wheels.

### 1.2 Non-goals for phase 1

These are on the roadmap in section 13 but are out of scope here:

- line mixing (Rosenkranz Y/G) and the AER line-file loader;
- secondary broadeners (H₂, He, CO₂, …), although the internal design must allow them (section 5.1);
- speed-dependent Voigt, Hartmann–Tran and other non-Voigt profiles, and χ-factor sub-Lorentzian wings;
- DIT-style fast algorithms for line lists of 10⁶ lines or more, and lookup-table generation;
- collision-induced absorption, HITRAN `.xsc` cross sections and continua;
- downloading spectroscopic data;
- integration into sasktran2, although the API must make it possible (section 11);
- explicit SIMD; a `nightly` feature using `portable_simd` may be added later if benchmarks justify it (section 10).

## 2. Naming

| Item | Name |
|---|---|
| Rust core crate (crates.io) | `absrbrs` |
| Python distribution (PyPI) | `absrbrs` |
| Python import | `absrbrs` |
| Python extension module | `absrbrs._absrbrs` |
| GitHub repository | `usask-arg/absrbrs` |
| Local folder | `absrbrs` |

`absrbrs` was unregistered on both crates.io and PyPI on 2026-10-14.
The GitHub URL returned 404, which does not rule out a private repository.

## 3. Conventions

- **Units in the core API.**
  - Wavenumber: cm⁻¹.
  - Temperature: K.
  - Pressure: Pa.
  - Volume mixing ratio: dimensionless.
  - Cross section: m² per molecule. This matches sasktran2. The Python layer can return cm² per molecule on request.
- **Line parameters** use HITRAN units and conventions:
  - T_ref = 296 K;
  - S_ref in cm⁻¹/(molecule cm⁻²), including natural isotopic abundance;
  - widths and shifts in cm⁻¹ atm⁻¹, with widths as half widths at half maximum.
- **Per-molecule cross sections.** Because S_ref includes isotopic abundance, σ is per molecule of the species with all isotopologues present at natural abundance.
- **Absorber partial pressure.** p_self = VMR · p and p_air = (1 − VMR) · p.
- **Wavenumber grids** passed to the core must be strictly increasing. They need not be uniform.

## 4. Architecture

### 4.1 Workspace layout

```
absrbrs/                  virtual Cargo workspace
  Cargo.toml              members = ["absrbrs", "absrbrs-py"]
  absrbrs/                core crate (published; no pyo3, no ndarray)
    src/lines/            LineList, isotopologue metadata
    src/io/               HITRAN .par and HAPI .data/.header readers
    src/partition/        PartitionFunction trait, TIPS-2021 tables
    src/faddeeva/         w(z) and w′(z), Fast and Accurate tiers
    src/state/            line-state layer: q and ∂q/∂θ
    src/profile/          profile layer: windowed accumulation
    src/engine/           orchestration, tiling, threading
    src/params/           parameter identifiers and maps
    src/error.rs
    benches/  examples/  tests/
  absrbrs-py/             pyo3 cdylib → absrbrs._absrbrs (publish = false)
  python/absrbrs/         Python package (maturin python-source)
  tests/                  pytest
  reference/              separate uv project: Arb reference generator and data
  benchmarks/             separate uv project: comparisons with sasktran2 and HAPI
  docs/                   Sphinx + MyST-NB; docs/specs and docs/plans excluded from the site
```

The repository copies the tooling and conventions of scatrs (`~/dev/scat-rs`).
See section 12 for the details.

### 4.2 Two-layer derivative strategy

Every cross section is a sum of per-line contributions, and each contribution depends on the parameters only through four shape quantities:

    q = (A, ν₀′, σ_D, y)

- A is the line amplitude.
- ν₀′ is the shifted line centre.
- σ_D is the Doppler 1/e half width.
- y = γ_L / σ_D is the Voigt ratio.

The computation therefore splits into two layers.

1. **The line-state layer (`state`)** maps (condition, line parameters) to q and, when requested, ∂q/∂θ.
   - This is small scalar code that runs once per (condition, line).
   - It is generic over a minimal `Real` trait. Production runs it in f64 with hand-derived closed-form derivatives. Tests run the same code with dual numbers to check those closed forms.
2. **The profile layer (`profile`)** evaluates the Faddeeva kernel over the line's window on the grid. It accumulates:
   - values;
   - forward derivatives, as linear combinations of a kernel triplet;
   - adjoint dot products.

The chain rule connects the layers:
- in forward mode, ∂σ/∂θ = Σ_lines (∂σ/∂q)(∂q/∂θ);
- in adjoint mode, ∂L/∂θ = Σ_lines (∫ g ∂σ/∂q dν)(∂q/∂θ).

This split keeps the expensive per-grid-point work independent of the number of parameters, apart from a few multiply-adds per parameter.
It also lets later phases swap components without touching the rest:
- non-Voigt profiles add components to q;
- secondary broadeners change only the line-state layer;
- a DIT engine replaces only the profile layer.

### 4.3 Core API rules

These follow scatrs:

- **Buffers.** Outputs are written into caller-supplied, contiguous, row-major slices with explicit dimensions, which are checked. There is no ndarray or pyo3 in the public core API.
- **Threading.** No global thread pool. Calls take `Option<&rayon::ThreadPool>`; `None` runs serially.
- **Determinism.** Results are bitwise identical for any thread count (section 7.3).
- **Errors.** There is one `thiserror` error enum, marked `#[non_exhaustive]`. No finite input may cause a panic.
- **Safety.** `#![forbid(unsafe_code)]` and `#![warn(missing_docs)]`.
- **Toolchain.** Stable Rust. Inner loops are written over branch-free index ranges so the compiler can auto-vectorize them.
- **Size limits.** `MAX_*` limits on output sizes are checked before allocating.
- **Dependencies** are kept minimal: thiserror and rayon at first. Each further dependency needs a dependency audit.

## 5. Line-state layer (phase 1 physics)

### 5.1 Inputs

- **A condition** is (T, p, VMR).
- **Required line fields:** ν₀, S_ref, E″, γ_air, γ_self, n_air, δ_air, molecule and isotopologue.
- **Optional line fields and their defaults:**

  | Field | Default when absent |
  |---|---|
  | n_self | n_air |
  | δ_self | δ_air |
  | δ′_air (temperature coefficient of the air shift, cm⁻¹ atm⁻¹ K⁻¹) | 0 |

  With all optional fields at their defaults, the physics is identical to sasktran2's `LineAbsorber` without line mixing. When a database supplies these fields, the physics improves with no change to the API.

**Broadeners are represented internally as a list.** Each broadener carries:
- its partial pressure;
- its width, temperature exponent, shift and shift temperature coefficient.

Phase 1 builds the list from air and self only. Secondary broadeners are a later change to the data model and API, not to the algorithms.

### 5.2 Formulas

Pressures in the width and shift formulas are in atm (p_Pa / 101325). c₂ is the second radiation constant, 1.4387769 cm K.

- **Line intensity:**

      S(T) = S_ref · Q(T_ref)/Q(T) · exp(−c₂E″/T)/exp(−c₂E″/T_ref) · (1 − exp(−c₂ν₀/T))/(1 − exp(−c₂ν₀/T_ref))

- **Doppler 1/e half width:**
  - σ_D = (ν₀/c) · √(2 k_B T / m), where m is the isotopologue mass;
  - it uses the unshifted ν₀, as HAPI and sasktran2 do.
- **Lorentz half width:**

      γ_L = (T_ref/T)^{n_air} γ_air p_air + (T_ref/T)^{n_self} γ_self p_self

- **Shifted centre:**

      ν₀′ = ν₀ + (δ_air + δ′_air (T − T_ref)) p_air + δ_self p_self

- **Shape quantities:**
  - A = 10⁻⁴ · S(T) / (√π σ_D), where the factor 10⁻⁴ converts cm² to m²;
  - y = γ_L / σ_D.
- **Contribution of one line:**

      σ_line(ν) = A · K(x, y),   x = (ν − ν₀′)/σ_D,   K = Re w(x + iy)

  - The contribution applies within |ν − ν₀′| ≤ W and is zero outside that window.
  - With pedestal subtraction enabled, the contribution is A · (K(x, y) − K(x_W, y)), where x_W = W/σ_D. This is the MT-CKD convention.
  - The default W is 25 cm⁻¹, and pedestal subtraction is off by default.

### 5.3 Derivatives of q

∂q/∂θ is computed for exactly the requested parameter classes.

- **State parameters:**
  - T, with p and VMR fixed;
  - p, with T and VMR fixed;
  - VMR, with T and p fixed. This derivative includes the self-broadening and self-shift terms, which sasktran2 omits.
- **Line parameters:** ν₀, S_ref, E″, γ_air, γ_self, n_air, n_self, δ_air, δ′_air, δ_self.
- **Partition function.** The temperature derivative includes −d ln Q/dT. This comes from the TIPS interpolant's own derivative (section 5.5), so values and derivatives are exactly consistent.
- **Optional fields.** A derivative with respect to an optional field that is absent is taken at its default value. For example, ∂/∂n_self at n_self = n_air is still defined.

The key closed forms, with u = c₂ν₀/T:

    ∂ln S/∂T  = c₂E″/T² − (u/T) e^{−u}/(1 − e^{−u}) − d ln Q/dT
    ∂ln σ_D/∂T = 1/(2T)
    ∂ln A/∂T  = ∂ln S/∂T − 1/(2T)
    ∂y/∂θ     = (∂γ_L/∂θ − y ∂σ_D/∂θ) / σ_D

### 5.4 Line selection and culling

- Lines are stored sorted by ν₀.
- For each grid tile, the candidate lines are found by binary search over [tile_min − W − Δ, tile_max + W + Δ]. Δ bounds the largest pressure shift over the conditions in the call, so shifted lines are never missed.
- Intensity culling is off by default. When it is enabled:
  - the cull decision is made per (condition, line) and treated as fixed for that call;
  - derivatives omit culled lines and add no terms for the cull threshold.

### 5.5 Partition functions

- **The `PartitionFunction` trait** returns Q(T) and d ln Q/dT for a (molecule, isotopologue) pair. It also returns `MissingPartitionFunction` or `OutsideTable` errors.
- **Built-in implementation.** TIPS-2021 tables, tabulated at 1 K steps.
  - Interpolation uses a smooth piecewise polynomial whose derivative is used directly, rather than a finite difference.
  - At table nodes, Q must equal HAPI's `partitionSum` exactly.
  - Between nodes, Q must agree with HAPI to within the interpolation error, which is measured and documented.
- **Table storage.** The tables are generated by a committed script that records provenance and a checksum. They are stored as a compact binary blob in the crate, with a size target of 3 MB or less; crates.io's limit is 10 MB.
- **Licence gate.** The redistribution terms of the TIPS-2021 data and HAPI must be checked before the tables are committed (license-compliance).
  - If redistribution is not allowed, the fallback has two parts: user-supplied tables, and a Python helper that builds the tables from an installed HAPI at runtime.
- **User-supplied tables.** Users can supply their own partition functions as per-isotopologue (T, Q) tables. Python callbacks are deliberately not supported, so the computation never needs the GIL.

### 5.6 Line data

- **`LineList`** stores structure-of-arrays columns for the required and optional fields, sorted by ν₀. It also holds isotopologue metadata: mass, plus abundance for reference.
- **Readers:**
  - HITRAN 160-character `.par` records;
  - HAPI `.data`/`.header` caches. The header gives the extra columns (n_self, δ_self and so on) when they are present. This lets users reuse the caches that sasktran2 already downloads.
- **Parsing rules:**
  - isotopologue codes 0, A and B map to 10, 11 and 12, as in sasktran2;
  - records with S = 0 are skipped;
  - a malformed record produces `Parse { path, line_number, message }` and never a panic.
- **Construction from arrays.** Lines can also be built directly from arrays (`LineList::from_columns`). A modified copy with a replaced column can be made cheaply, which is how an optimizer perturbs line parameters.

## 6. Faddeeva kernel and profile layer

### 6.1 Kernel contract

- Every evaluation returns w(z) and w′(z) together, using the same formula in each region.
  - The derivative is never recovered from the identity w′ = −2zw + 2i/√π in regions where that identity cancels catastrophically. In the wings this cancellation loses about |z|² in relative accuracy.
  - Using one formula for both means value and derivative stay mutually consistent, which matters for optimizers.
- **The kernel triplet** used by the profile layer is:
  - K = Re w;
  - K_x = ∂K/∂x = Re w′;
  - K_y = ∂K/∂y = −Im w′, by the Cauchy–Riemann equations.

### 6.2 Validated domain and error metric

- **Domain:**
  - y ∈ [0, 10⁶]. y = 0 is zero pressure, a pure Gaussian. The reference grid includes y = 0 exactly and is log-spaced from 10⁻³⁰ upwards. Tiny nonzero y needs testing in its own right: far from the centre the Lorentz tail y/(√π x²) dominates exp(−x²) even when y is tiny;
  - |x| ≤ 10⁸. This covers 25 cm⁻¹ windows on microwave lines, where σ_D is about 10⁻⁶ cm⁻¹.
- **Value error** is measured relative to K, which is positive everywhere. Points where the true K is below the smallest normal f64 are exempt.
- **Derivative error** is measured for each component (K_x and K_y), relative to |w′(z)|. This avoids dividing by zero where K_x crosses zero at the line centre.

### 6.3 Tiers

| Tier | Value error | Derivative error | Algorithm |
|---|---|---|---|
| **Fast** (default) | ≤ 1 × 10⁻⁶ | ≤ 1 × 10⁻⁵ | Chosen by measurement in phase 0. The candidates are combined region by region: a Weideman N = 24/32 rational core, a Humlíček region-I asymptotic for large \|z\|, a higher-order asymptotic tail, and a small-y expansion for the hard region at small y and 3 < \|x\| < 30. |
| **Accurate** | ≤ 2 × 10⁻¹³ | ≤ 2 × 10⁻¹² | A Rust port of Steven G. Johnson's MIT-licensed Faddeeva Package: continued fractions for large \|z\|, and erfcx/Dawson-based evaluation elsewhere. ACM TOMS Algorithm 916 code is not used because of the ACM software licence. |

- **Choosing the Fast algorithm.** Phase 0 builds a Python prototype of the candidates and measures their accuracy maps against Arb. The Fast-tier algorithm is the cheapest combination that meets both targets over the whole domain.
- **Reference point.** sasktran2's 8-pole J-pole approximation does not meet the Fast target. Its relative error a few Doppler widths from the line centre is about 10⁻³ at y = 10⁻³ and above 100% at y = 10⁻⁶. See Appendix A.

### 6.4 Profile accumulation

For each line, the region thresholds in x are converted to grid-index ranges by binary search. Each range is then evaluated with a branch-free loop.

- **Value:** σ[i] += A · K(xᵢ, y).
- **Forward derivative** for each requested parameter θ:

      ∂σ[i]/∂θ += a_θ K + (b_θ + c_θ xᵢ) K_x + d_θ K_y + e_θ

  The per-line coefficients are:
  - a_θ = ∂A/∂θ;
  - b_θ = −A (∂ν₀′/∂θ)/σ_D;
  - c_θ = −A (∂σ_D/∂θ)/σ_D;
  - d_θ = A ∂y/∂θ;
  - e_θ is the pedestal term, which is constant over the window and zero when pedestal subtraction is off.

  This costs about four fused multiply-adds per parameter per grid point, on top of one kernel evaluation.
- **Specialized state-only kernels.** The state-only path is monomorphized for 0, 1, 2 or 3 state directions, so it never computes line-parameter coefficients.
- **Adjoint.** For a weight vector g (∂L/∂σ), each line accumulates over its window:
  - G₀ = Σ g K;
  - G₁ = Σ g K_x;
  - G₂ = Σ g x K_x;
  - G₃ = Σ g K_y;
  - G₄ = Σ g, only when pedestal subtraction is on.

  These give ∂L/∂q, and ∂L/∂θ then follows for every θ through ∂q/∂θ.

## 7. Engine and Rust API

### 7.1 Shape

```rust
let lines = absrbrs::io::read_hitran_par("CO2.par")?;
let pf    = absrbrs::partition::Tips2021::builtin();
let opts  = Options::default()
    .with_tier(Tier::Fast)
    .with_cutoff(Cutoff { width_cminv: 25.0, pedestal: false });
let lbl   = LineByLine::new(&lines, &pf, opts)?;

lbl.cross_section(&grid, &conditions, out, pool)?;                          // out: [n_cond × n_grid]
lbl.cross_section_jacobian(&grid, &conditions, &params, out, jac, pool)?;   // forward
lbl.cross_section_vjp(&grid, &conditions, grad, &mut line_grad, pool)?;     // adjoint
```

- `conditions` is `&[Condition { temperature_k, pressure_pa, vmr }]`.
- `LineByLine::new` validates the line list, resolves per-isotopologue masses and checks that a partition function exists for every isotopologue present.
- Exact names are settled in the implementation plan. The structure above is fixed.

### 7.2 Parameters (`params`)

- **State parameters** (T, p, VMR) act per condition.
  - Their Jacobian block has shape [n_state × n_cond × n_grid].
  - There are no cross-condition terms.
- **Line parameters** are linear combinations of per-line fields, ∂σ/∂θⱼ = Σ_(l,f) w_(l,f,j) ∂σ/∂p_(l,f).
  - Their Jacobian block has shape [n_line_params × n_cond × n_grid].
  - Convenience constructors cover the common cases:
    - a single field of a single line, with weight 1;
    - a multiplicative scale factor on a field over a set of lines, evaluated at s = 1, with weights equal to the field values;
    - an additive offset on a field over a set of lines.
- **Adjoint output.** `LineGradient` holds:
  - one length-n_lines array per line field, summed over conditions;
  - an [n_cond × 3] array of state gradients.

### 7.3 Work decomposition and determinism

- **Line state** (q and, when needed, ∂q/∂θ) is computed once per (condition, line). It is held only for the conditions currently in flight.
- **Tiles.** The grid is divided into tiles of a fixed size. The size is a constant and does not depend on the thread count. Each task is one (condition, tile) pair on the caller's pool.
- **Summation order.** Each grid point sums its lines in sorted-ν₀ order, so values and forward derivatives are bitwise independent of thread count and tiling.
- **Adjoint sums.** Each line's adjoint partial sums are combined across tiles in tile order. Because tile boundaries are fixed, the adjoint is also bitwise deterministic across thread counts.

### 7.4 Errors

| Variant | Raised when |
|---|---|
| `InvalidInput(String)` | the grid is not strictly increasing; T ≤ 0; p < 0; VMR is outside [0, 1]; an input is NaN or infinite; a buffer shape is wrong; an output size exceeds its limit |
| `MissingPartitionFunction { molecule, isotopologue }` | an isotopologue has no partition function |
| `OutsideTable { molecule, isotopologue, temperature }` | T is outside the partition-function table |
| `Parse { path, line_number, message }` | a line-file record is malformed |
| `Io { path, kind, message }` | a file cannot be read; `kind` is the `std::io::ErrorKind`, so a missing file can be told from other failures |

## 8. Validation

### 8.1 Arb reference suite (`reference/`)

- `reference/` is a separate uv project using python-flint, as in scatrs.
- It evaluates w(z) = exp(−z²) erfc(−iz) and w′(z) in arbitrary precision on a log-spaced (x, y) grid. The grid includes y = 0 and the boundaries of the domain in section 6.2.
- Values are committed as JSONL, together with a manifest and error bounds. CI regenerates the data and fails if anything differs.
- The suite also produces brute-force arbitrary-precision cross sections, with derivatives, for small synthetic line lists. These cover several conditions, pedestal on and off, and windows that straddle the edges of the grid.

### 8.2 Accuracy maps

- Accuracy maps of both tiers, for values and for derivatives, are computed against the Arb grid.
- They are enforced in CI and published in the docs.

### 8.3 Invariant tests (proptest)

- K(−x, y) = K(x, y), and the corresponding parities of K_x and K_y.
- K > 0, and K is non-increasing in |x|.
- ∫K dx = √π, checked numerically.
- Cross sections scale linearly with S_ref.
- Results are bitwise identical at 1, 2 and 8 threads.
- The parsers never panic on arbitrary input.

### 8.4 Derivative verification

- **Line-state closed forms vs dual numbers:** agreement to about 1 × 10⁻¹³.
- **Forward vs finite differences:** cross-section derivatives agree with Richardson-extrapolated finite differences, for every state and line parameter.
- **Forward vs adjoint:** they are checked against each other with a dot-product test, ⟨g, J v⟩ = ⟨Jᵀ g, v⟩, to about 1 × 10⁻¹².
- **Specialized vs general path:** the state-only path agrees with the general path to within 1 × 10⁻¹⁵.
- **Partition function:** d ln Q/dT is checked against finite differences of the interpolant.

### 8.5 Cross-code comparisons (`benchmarks/`)

These run in a scheduled job with pinned versions, and the results are published in the docs.

- **HAPI** `absorptionCoefficient_Voigt`, with matched cutoffs, partition functions and line selection.
- **sasktran2** `LineAbsorber`, with the sasktran2-matching defaults.
  - The expected differences are documented: J-pole accuracy, the Gaussian-region test (Appendix A), and the VMR Jacobian.
- **Speed comparison** against sasktran2 on real HITRAN bands, covering values and T/p/VMR Jacobians.

### 8.6 Test data

- Small HITRAN snippets are used in the loader golden tests.
- HITRAN's redistribution terms must be checked before any snippet is committed.
- Rust unit tests and benchmarks otherwise use deterministic synthetic line lists.

## 9. Python layer

### 9.1 API

```python
import absrbrs as ab

lines = ab.LineList.from_hitran("CO2.par")              # or ab.LineList.from_hapi(cache_dir, "CO2")
lines = lines.select(isotopologues=[1, 2], wavenumber=(6200.0, 6400.0))

xs = ab.cross_section(
    lines, wavenumber=grid, temperature=T, pressure=p, vmr=vmr,
    tier="fast", cutoff=25.0, pedestal=False,
    derivatives=["temperature", "pressure", "vmr"],
    num_threads=4,
)
xs.values            # broadcast(T, p, vmr).shape + (n_grid,)
xs.d("temperature")
xs.to_xarray()       # requires the optional xarray extra
```

- **`ab.LineByLine(lines, tier=..., cutoff=..., pedestal=...)`** is the reusable calculator for repeated calls, for example inside an optimizer.
  - It provides `.cross_section(...)` and `.vjp(...)`.
  - `ab.cross_section` is a thin wrapper around it.
- **Line-parameter derivatives (opt-in):**
  - `line_derivatives=[ab.line_param("S", lines=[10, 11]), ab.line_scale("gamma_air", lines=mask), ...]`;
  - each one is read back with `xs.d(name)`.
- **Adjoint:** `lbl.vjp(grad, wavenumber=..., temperature=..., pressure=..., vmr=...)` returns a `LineGradient` with:
  - one array per line field (`.S`, `.gamma_air`, …);
  - condition-shaped `.temperature`, `.pressure` and `.vmr`.
- **`LineList`** is a pyclass that owns the Rust columns, so nothing is copied per call.
  - Columns are exposed as read-only numpy arrays.
  - `LineList.from_arrays(...)` builds a list from arrays.
  - `lines.with_column(name, values)` returns a modified copy.
- **Partition functions:**
  - `ab.partition.tips2021(molecule, isotopologue, T)` returns Q and d ln Q/dT;
  - custom partition functions are passed as `{(molecule, isotopologue): (T, Q)}` tables.
- **Unsorted grids** are sorted in Python and the results are returned in the caller's original order, as sasktran2 does.

### 9.2 Binding conventions (from scatrs)

- Python converts types and shapes; Rust validates values.
- The GIL is released for the whole computation.
- Parameters use `Bound<PyArray…>` with `try_readonly()`; a failed borrow raises `BufferError`. `PyReadonlyArray` is never used as a parameter type.
- Outputs are allocated by the binding, after a size check.
- **Exceptions:**
  - `AbsrbrsError` is the base class;
  - `InvalidInputError` subclasses both `AbsrbrsError` and `ValueError`;
  - `ParseError` and `OutsideTableError` subclass both `AbsrbrsError` and `ValueError`;
  - `MissingPartitionFunctionError` subclasses both `AbsrbrsError` and `LookupError`;
  - a file that cannot be read raises the built-in `FileNotFoundError`, `PermissionError` or `OSError`, chosen by the `std::io::ErrorKind` of `Error::Io`.
- The extension has hand-written `.pyi` stubs and ships `py.typed`. CI checks them with `mypy.stubtest`.

### 9.3 Packaging

- pyo3 and numpy 0.29, matching scatrs.
- abi3 wheels for Python 3.11 and later.
- numpy ≥ 1.26.
- An optional `xarray` extra.
- The version comes only from the Cargo workspace (`dynamic = ["version"]`).

## 10. Performance

### 10.1 Benchmark workload

The primary workload is atmospheric radiative transfer:
- 100 conditions sharing one grid;
- 10⁵ grid points;
- synthetic line lists sized like the O₂ A-band (about 10³ lines) and the CO₂ 1.6 µm bands (about 3 × 10⁴ lines);
- values only, and values plus T/p/VMR Jacobians.

Criterion benchmarks run in Rust. Comparisons against sasktran2 on real data live in `benchmarks/`.

### 10.2 Targets

| Quantity | Target |
|---|---|
| Fast-tier values | At least as fast as sasktran2's value path, at the same thread count and cutoff |
| Values + 3 state derivatives | ≤ 1.5 × our own value path |
| Accurate tier | ≤ 4 × Fast tier |
| Parallel efficiency | ≥ 80% at 8 threads on the primary workload |

- **The sasktran2 parity target is the main risk.** sasktran2 uses nightly f64x8 SIMD. If auto-vectorized stable code misses parity, the remedy is an optional `nightly` cargo feature using `portable_simd` in the kernel loops, evaluated at that point.
- **Line-parameter workloads** are reported, but their speed does not gate any phase.

### 10.3 Build profile

- Release builds use `lto = "fat"` and `codegen-units = 1`.
- The test profile uses `opt-level = 2`.

## 11. Constraints for future sasktran2 integration

- **No pyo3 or ndarray in the core public API.** sasktran2 uses pyo3 0.27 and ndarray 0.16; this crate uses pyo3 0.29 only in `absrbrs-py`.
- **Threading.** The caller supplies the rayon pool. sasktran2 owns its pool through `Config.num_threads`.
- **Layout and units match sasktran2:** σ as [n_cond × n_grid] in m², wavenumber in cm⁻¹, and Jacobian keys corresponding to `temperature_k` and `pressure_pa`.
- **The VMR derivative is new.** sasktran2's VMR weighting function currently ignores self-broadening. When integrating, the constituent should use the new VMR derivative.
- **sasktran2 options map directly:**
  - `line_contribution_width` maps to `Cutoff::width_cminv`;
  - `cull_factor` maps to the opt-in culling.
  - `line_coupling` waits for phase 2.
- **Toolchain.** The stable-Rust code builds unchanged under sasktran2's pinned nightly toolchain.

## 12. Build, CI, release, docs

The setup copies scatrs as of 2026-10-14. The workflows are:

- **`ci.yml`:**
  - lint: `cargo fmt`, `clippy -D warnings`, rustdoc `-D warnings`, `cargo package` followed by testing the packaged crate, ruff and actionlint;
  - Rust tests on Ubuntu, macOS and Windows;
  - an MSRV check (Rust 1.85, edition 2024);
  - Python wheel tests, including a numpy 1.26 job;
  - a strict Sphinx build plus `mypy.stubtest`;
  - a check that regenerated reference data matches the committed data;
  - a single `all-checks` job as the required check.
- **`release.yml`:** abi3 wheels for linux, musllinux, windows and macOS, plus an sdist. Publishing to PyPI and crates.io uses Trusted Publishing on `v*` tags.
- **`docs.yml`:** GitHub Pages.
- **`scheduled.yml`:** reference drift, Intel macOS, ignored (slow) tests, and the cross-code comparisons.

Supporting conventions:

- Every action is pinned to a commit SHA.
- `persist-credentials: false` on checkout.
- dependabot with grouped updates.
- uv for Python tooling, with `uv.lock` committed.
- pre-commit hooks.
- `CHANGELOG.md` in Keep a Changelog format.
- MIT licence.
- `CONTRIBUTING.md` lists the exact commands CI runs.

**Docs** use Sphinx with furo and MyST-NB, executing notebook examples. The sections are:
- user guide;
- theory: line shapes, Faddeeva tiers and the Jacobian mechanics;
- validation: accuracy maps and the HAPI and sasktran2 comparisons;
- API reference.

## 13. Build order and roadmap

| Phase | Contents |
|---|---|
| **0** | Repository, CI, Python package skeleton and reference-generator scaffolding. A Python Faddeeva prototype and Arb accuracy maps to choose the Fast-tier algorithm. Licence checks for TIPS-2021, HAPI, HITRAN snippets and the Faddeeva Package. |
| **1a** | Faddeeva kernels at both tiers, with derivatives and accuracy maps in CI. |
| **1b** | `LineList`, the HITRAN and HAPI readers, TIPS-2021 partition functions, and the line-state layer with closed-form ∂q/∂θ checked by dual numbers. |
| **1c** | Profile layer and engine: values, the state-only fast path, forward line parameters, adjoint, tiling and determinism, plus the Python API. |
| **1d** | Cross-code validation against HAPI and sasktran2, benchmarks against the performance targets, docs, and release 0.1.0. |
| 2 | Line mixing (Rosenkranz Y/G with temperature interpolation), the AER v3.8.1 loader, secondary broadeners. |
| 3 | Speed-dependent Voigt and Hartmann–Tran profiles; χ-factor and sub-Lorentzian wings. |
| 4 | DIT-style algorithms for very large line lists (HITEMP, ExoMol); lookup-table generation. |
| 5 | Collision-induced absorption, HITRAN `.xsc` cross sections, continua. |
| 6 | sasktran2 integration: an adapter in sasktran2 replacing `LineAbsorber`. |
| 7 | Data access: downloading from HITRAN, HITEMP and ExoMol. |

Each phase follows the scatrs process: a dated spec when needed, then a dated plan, then a follow-ups file.

## Appendix A. Observations about the current sasktran2 line-by-line code

These were found while reviewing `~/dev/sasktran2` at commit 70d6bb8a. The first six were checked directly in the source or by measurement. The seventh was reported during code exploration and still needs checking.

1. **The Gaussian-only region test compares raw wavenumbers.** In `split_and_assign` (`rust/sasktran2-rs/src/optical/types/line_absorber.rs:163`), the test compares max(|ν_min|, |ν_max|), in absolute cm⁻¹, against a threshold in units of x = (ν − ν₀)/σ_D. As a result the pure-Gaussian branch essentially never runs. The derivative path mirrors the same test.
2. **J-pole accuracy.** The 8-pole approximation in `math/errorfunctions/optimized.rs` was compared with `scipy.special.wofz` on 2026-10-14. Its maximum absolute error is about 10⁻⁷ of the peak. Its maximum relative error, near |x| ≈ 10, is:

   | y | max relative error |
   |---|---|
   | 10⁻⁶ | 1.3 |
   | 10⁻³ | 1.3 × 10⁻³ |
   | 10⁻² | 1.3 × 10⁻⁴ |
   | 1 | 9 × 10⁻⁷ |

3. **Self broadening uses n_air.** The same temperature exponent n_air is applied to both air and self broadening (`optical/line/db.rs`, `adjusted_parameters`).
4. **The pressure shift is δ_air · p_total.** It has no temperature dependence and no separate self shift.
5. **`subtract_pedestal` is never used.** It is stored and has a builder method, but the computation never reads it.
6. **`between_slice` scans linearly.** It is O(N_lines) instead of a binary search (`optical/line/db.rs:160`).
7. **The VMR weighting function omits self broadening** (reported; to be checked). It is reported to be computed as σ · N, without the term from the dependence of σ on VMR through self broadening.
