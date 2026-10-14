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
    7th significant figure plus the f32 rounding instead.
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

## Measured

- **Readers** (the opt-in `ABSRBRS_HITRAN_LINES_DIR` test of Task 3, on the user's local HAPI
  caches): every real cache parses in at most 0.1 s; the largest, CO2, has 545 084 lines.
  Through Python, `LineList.from_hapi` on CO2 takes 0.11 s, against 0.10 to 0.12 s in Rust; the
  GIL is released while reading and while evaluating TIPS.
- **TIPS package size.** `absrbrs/data/tips2021.bin` is 3 508 680 bytes (181 isotopologues,
  876 441 f32 values); the packaged crate (`cargo package -p absrbrs`) is 3 321 050 bytes,
  about 3.3 MB.
- **Line state against Arb** (`reference/data/line_state.jsonl`, 192 records at 256 bits): the
  largest relative error of A, ν₀′, σ_D and y is 4.0e-15 (in A), against a tolerance of 1e-13.
  On 360 cases the closed-form Jacobians match dual numbers to within 0.057 of the bound
  |closed − dual| ≤ 1e-13 (|closed| + |dual|).
- **Line-state speed** (release, M-series Mac): `line_shape` 21 ns, the state Jacobian 27.5 ns
  and the line Jacobian 42 ns per line and condition.
- **TIPS-2021 against HAPI** (`validation/tests/test_tips.py`, hitran-api 1.3.0.0,
  `version=2021`): 144 isotopologues, every HAPI node from 10 K to 1000 K, 14 396 nodes in all.
  HAPI 1.3.0.0's TIPS-2021 has no table for 56_1 or 57_1, which absrbrs has; the other 10 HAPI
  isotopologues without a comparison have no TIPS-2021 table in either. The 10 K and 20 K nodes
  of 31_2 and 31_3 are skipped (4 nodes), being below the start of absrbrs's tables; no other
  table loses a node. Every HAPI value carries at most 7 significant figures. The largest
  relative difference is 5.5e-7 (4_3 at 560 K: absrbrs 10018.1445, HAPI 10018.15), and every
  isotopologue has some node beyond 2e-7. Against the bound of half a unit in HAPI's 7th
  significant figure plus the f32 rounding, the worst node is at 0.992 of it (20_2 at 550 K).

## TIPS-2021 source defects

Found while building the tables (`tools/tips2021.py`); worth reporting upstream.

- **H₂S 31_2 and 31_3 are unphysical below about 30 K.** Their ratios to 31_1 should be close
  to the ratios of the nuclear-spin degeneracies (about 1 and 4), but are 0.44 and 1.77 at
  2 K, 1.43 and 5.72 at 6 K and 1.03 and 4.11 at 20 K; from 30 K they are within 0.6%.
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
  (`value()` costs about 1.3 ns per call). The engine resolves masses with
  `lines::isotopologue` and partition values with `PartitionFunction`, and computes
  `PartitionRatio` once per (condition, isotopologue).
- Hoist the per-line constants into that precompute. The line-state functions recompute, per
  line and condition, quantities that depend on the line only: expm1(−c₂ν₀/T_ref) and
  (ν₀/c)√(2k_B/m). `line_shape_line_jacobian` also recomputes the stimulated-emission factor,
  the Boltzmann factor and the intensity factor that `line_shape` has already computed (for
  ∂A/∂S_ref, in `absrbrs/src/state/jacobian.rs`).
- The engine checks the line-state preconditions (T > 0, mass > 0, 0 ≤ VMR ≤ 1); the layer
  itself does not.

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
