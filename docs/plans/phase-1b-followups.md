# Phase 1b follow-ups and decisions

## Decisions

- **TIPS-2021 at 1 K steps, as f32.** The built-in tables are the 1 K `.QTpy` tables of the
  Zenodo record (doi:10.5281/zenodo.4708099), stored as f32 (≤ 6e-8 relative rounding). This
  supersedes spec §5.5's "equal HAPI's partitionSum exactly at nodes". The checks are now:
  - against the source, to the f32 rounding (`absrbrs/src/partition/mod.rs`), and against the
    Fortran tables of `BD_TIPS_2021.FOR` to their 5 significant figures at all 440 011 shared
    nodes (worst 0.9988 of the bound, recorded in `absrbrs/data/tips2021.json`);
  - against HAPI at its 10 K nodes, to HAPI's rounding (see "Measured"). The plan's 2e-7
    relative tolerance could not hold, because HAPI rounds Q to 7 significant figures, which
    alone moves it by up to 5e-7 relative; the test bounds each node by half a unit in HAPI's
    7th significant figure plus the f32 rounding instead;
  - against HAPI between its nodes (spec §5.5), at the non-node integer temperatures from
    150 K, to 1.2e-5 (see "Measured"). This measures HAPI's 10 K interpolation against the
    source, and is why HAPI's tables are not the reference.

  37 of the 181 tables (the 35 whose isotopologues are not in `hapi.ISO`, plus 56_1 and 57_1,
  which HAPI's TIPS-2021 lacks) are checked only against the Fortran tables.
- **Interpolation.** TIPS-2021 is interpolated in ln Q against T, user tables
  (`TabulatedPartition`) in ln Q against ln T (exact for Q ∝ Tⁿ). Both use cubic Hermite
  polynomials with numpy-`gradient`-style second-order node slopes (`edge_order=2`; the secant
  for two-point tables), and d ln Q/dT is the interpolant's own derivative.
- **Trimmed TIPS tables.** H₂S 31_2 and 31_3 start at 30 K and CH₃ 57_1 at 2 K (see "TIPS-2021
  source defects"); a temperature below a table's start is `OutsideTable`.
- **The tie rule for absent optional fields** refines spec §5.3: the derivative with respect to
  an absent field is the partial at its default, and the derivative with respect to the field
  it defaults to includes its effect through the default (∂/∂n_air contains the self term when
  n_self is absent, ∂/∂δ_air the self shift when δ_self is absent). This keeps the Jacobian
  consistent with finite differences of `with_column`.
- **The constants are exact SI**, with c₂ = hc/k_B = 1.4387768775039338 cm K, the nearest
  double; sasktran2's 1.4387769 is 1.6e-8 relative larger.
- **The line-state layer stays crate-private** until the phase 1c engine uses it.
- **Isotopologue masses** come from one global table (`lines::isotopologue`), not from metadata
  held in each `LineList` as spec §5.6 has it; custom isotopologues and overridden masses are
  not possible yet.
- **Python column accessors return read-only copies**, not zero-copy views (spec §9.1 says
  nothing is copied per call). That is acceptable at line-list sizes; revisit if profiling shows
  otherwise.
- **`with_column("nu", …)` re-sorts the list** (documented); per-line arrays in the old order
  no longer line up.
- **Readers.** Both skip records with ν₀ ≤ 0, which a `LineList` cannot hold, as well as
  records with S_ref = 0 (as sasktran2 does). Negative S_ref and negative widths are `Parse`
  errors carrying the line number, not `InvalidInput` from `from_columns`. Numbers may take
  the Fortran forms that HITRAN files use; a stray CR (other than before the LF) is rejected.
  An optional field listed in both `order` and `extra` of a HAPI header is read from its
  fixed-width column (`order` is preferred).
- **`LineField` is `#[non_exhaustive]`**, with `LineField::COUNT` (10) for arrays indexed by it,
  and `LineList` equality treats NaN optional entries as equal, so that a list equals its clone.
- **`tips2021()`** (a function returning the static tables) replaces spec §7.1's
  `Tips2021::builtin()`.
- **Broadeners.** Spec §5.1's broadener list is not used yet: air and self are hardcoded in the
  crate-private `LineParameters`.
- **serde_json is a core dependency** (it parses HAPI's JSON headers); it brings
  serde_core, itoa, memchr and zmij, all MIT or dual-licensed. They are added to the wheel-notice
  follow-up in `docs/plans/phase-1a-followups.md` ("Before the first release").
- **`LineList::value` and `has_value` panic on an out-of-range index**, against spec §4.3 (no
  panics on finite input). Phase 1c decides between a documented carve-out and `Option`
  accessors.
- **`_absrbrs.isotopologue_info`** has no public Python wrapper yet; it is reserved for phase 1c.

## Measured

- **Readers** (the opt-in `ABSRBRS_HITRAN_LINES_DIR` test of Task 3, on the user's local HAPI
  caches): every real cache parses in at most 0.1 s; the largest, CO2, has 545 084 lines.
  Through Python, `LineList.from_hapi` on CO2 takes 0.11 s, against 0.10 to 0.12 s in Rust; the
  GIL is released while reading and while evaluating TIPS.
- **TIPS package size.** `absrbrs/data/tips2021.bin` is 3 508 680 bytes (181 isotopologues,
  876 441 f32 values); the packaged crate (`cargo package -p absrbrs`) is about 3.3 MB.
- **Line state against Arb** (`reference/data/line_state.jsonl`, 192 records at 256 bits): the
  largest relative error of A, ν₀′, σ_D and y is 4.0e-15 (in A), against a tolerance of 1e-13.
  On 360 cases the closed-form Jacobians match dual numbers to within 0.057 of the bound
  |closed − dual| ≤ 1e-13 (|closed| + |dual|).
- **Line-state speed** (ad hoc, release build, Apple M-series; a criterion bench in phase 1c):
  `line_shape` 21 ns, the state Jacobian 27.5 ns and the line Jacobian 42 ns per line and
  condition.
- **TIPS-2021 against HAPI** (`validation/tests/test_tips.py`, hitran-api 1.3.0.0,
  `version=2021`): 144 isotopologues, every HAPI node from 10 K to 1000 K, 14 396 nodes in all.
  HAPI 1.3.0.0's TIPS-2021 has no table for 56_1 or 57_1, which absrbrs has; the other 10 HAPI
  isotopologues without a comparison have no TIPS-2021 table in either. The 10 K and 20 K nodes
  of 31_2 and 31_3 are skipped (4 nodes), being below the start of absrbrs's tables; no other
  table loses a node. Every HAPI value carries at most 7 significant figures. The largest
  relative difference is 5.5e-7 (4_3 at 560 K: absrbrs 10018.1445, HAPI 10018.15), and every
  isotopologue has some node beyond 2e-7. Against the bound of half a unit in HAPI's 7th
  significant figure plus the f32 rounding, the worst node is at 0.992 of it (20_2 at 550 K).
- **TIPS-2021 against HAPI between its nodes** (same test): every non-node integer temperature
  from 150 K to the end of HAPI's table, 624 969 temperatures. Away from the kinks the worst
  is 1.115e-5 (SF₆ 30_1 at 265 K), then 7.9e-6 (20_1 at 632 K) and 5.1e-6 (35_2 at 155 K); the
  bound is 1.2e-5. Within 20 K of the kinks HAPI's stencil spans the slope discontinuity:
  2.0e-3 (SO₃ 47_1 at 651 K), 2.2e-4 (23_2 at 1001 K) and 4.7e-5 (23_1 at 996 K), bounded by
  2.5e-3. Below 150 K HAPI's error grows, to about 9e-5 near 50 K (45_2 at 55 K).
- **Interpolation figures** (`uv run --locked --directory validation python -m
  absrbrs_validation.tips_interpolation`, quoted in `docs/theory/partition-functions.md` and
  the `Tips2021` rustdoc):
  - HAPI (3- and 4-point Lagrange in Q on 10 K) against the 1 K source at non-node integer T:
    at most 1.1e-5 at 200 to 300 K (SF₆ 30_1); under 1e-6 at 150 to 350 K for molecules 1 to 7,
    isotopologue 1 (worst O₃, 9.5e-7); up to 1.03 below 20 K (NO⁺ 36_1).
  - absrbrs against `TIPS_2021.py`'s linear interpolation in Q (1 K grid), at half-integer T:
    median 4.1e-6 and worst 2.7e-5 (30_1) at 150 to 350 K; worst 6.4e-5 from 50 K (35_1);
    1.4e-4 at the SO₃ kink.
  - absrbrs against 6-point Lagrange interpolation of ln Q, at half-integer T: median 2.7e-8
    and worst 3.6e-7 from 50 K, more than 4 K from the kinks; 5.1e-5 near the SO₃ kink. Below
    20 K, up to 1.05e-2 (36_1), and up to 0.16 against linear in Q (0.58 for 57_1): the
    "1e-2 to 2e-1 below 20 K" of the docs is this spread between interpolants, not a measured
    error.

## TIPS-2021 source defects

Found while building the tables (`tools/tips2021.py`); worth reporting upstream.

- **H₂S 31_2 and 31_3 are unphysical below about 30 K.** Their ratios to 31_1 should be close
  to the ratios of the nuclear-spin degeneracies (about 1 and 4), but are 0.44 and 1.77 at
  2 K, 1.43 and 5.72 at 6 K and 1.03 and 4.11 at 20 K. From 30 K to about 300 K they are
  within 0.6% (at most 0.58% and 0.46%); above that they drift, as expected physically, to
  1.5% at 1000 K and 5.5% at 4000 K.
- **Q(1 K) ≤ 0 for 31_2, 31_3 and 57_1** (negative for the two H₂S tables, 0 for CH₃), where
  ln Q is undefined. 57_1's low-temperature values also carry few digits (Q(2 K) = 0.000416).
- **CF₄ 42_1's maximum.** `TIPS_2021.py` gives Tmax = 3010 K, but the table runs, smoothly, to
  5000 K; absrbrs uses all of it. (HAPI's TIPS-2021 table for 42_1 stops at 3010 K.)

## Review notes

- **Python `LineList`:** no pickling (so no multiprocessing) and no `__eq__`. Every property
  access copies; if loops over `lines.nu[i]` become common, cache the read-only arrays per
  instance.
- **Python `partition.tips2021`** reads a view of the caller's array with the GIL released, the
  same accepted pattern as in phase 1a.
- **`OSError`'s `errno` and `filename` are `None`**, by the phase 1a `to_pyerr` design.

## For phase 1c

- Build each line's condition-independent `LineParameters` once in `LineByLine::new`
  (`value()` costs about 1.3 ns per call: ad hoc, release build, Apple M-series; a criterion
  bench in phase 1c). The engine resolves masses with
  `lines::isotopologue` and partition values with `PartitionFunction`, and computes
  `PartitionRatio` once per (condition, isotopologue).
- Hoist the per-line constants into that precompute. The line-state functions recompute, per
  line and condition, quantities that depend on the line only: expm1(−c₂ν₀/T_ref) and
  (ν₀/c)√(2k_B/m). `line_shape_line_jacobian` also recomputes the stimulated-emission factor,
  the Boltzmann factor and the intensity factor that `line_shape` has already computed (for
  ∂A/∂S_ref, in `absrbrs/src/state/jacobian.rs`).
- The engine checks the line-state preconditions (T > 0, mass > 0, 0 ≤ VMR ≤ 1); the layer
  itself does not.
- 35 TIPS-2021 tables have no isotopologue mass, because `hapi.ISO` lacks them (1_8, 1_9,
  2_13, 3_6 to 3_18, 5_7 to 5_9, 7_4 to 7_6, 22_3, 27_3, 38_3, 41_2 to 41_4, 44_2 to 44_6,
  48_2, 54_2). Decide which error the engine raises for lines of these isotopologues.
- The engine needs Q(296 K) for the ratio, so a user table that does not cover 296 K would
  fail with a confusing `OutsideTable` at 296 K. Validate this in `LineByLine::new`, with a
  clear error.
- Python user partition tables (spec §9.1): `TabulatedPartition` is still Rust-only.
- The crate-private `state::Condition` will sit beside the public engine `Condition`; rename
  or merge them.

## For phase 1d (sasktran2 comparison)

HAPI 1.3's reader for tables with extras (`getRowObjectFromString`, used by `storage2cache`)
parses `local_iso_id` with `int()`. So `.par` code `0` becomes isotopologue 0 instead of 10,
and lines with codes `A` and `B` are silently dropped. sasktran2's `CO2.nc`, converted through
HAPI, therefore has 539 603 lines where the data has 545 084 (3 980 `A` and 1 501 `B` lines
missing), and 2 926 lines labelled isotopologue 0. absrbrs follows the HITRAN numbering, so a
comparison must expect the extra lines and isotopologues 10, 11 and 12. (sasktran2's own Rust
`.par` loader maps the codes correctly.)

## For phase 2

- The AER line lists (a local copy is in sasktran2's database).
- Line mixing.
