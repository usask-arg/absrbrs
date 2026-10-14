# Partition functions

The line intensity at temperature T needs the total internal partition sum Q(T) of the line's
isotopologue, and the temperature Jacobian needs d ln Q/dT. A partition function in absrbrs
returns both, from the same interpolant, so they are consistent with each other. A source of
partition functions implements the `PartitionFunction` trait; absrbrs has two, the built-in
TIPS-2021 tables and user tables.

## TIPS-2021

The built-in tables are TIPS-2021 (Gamache et al. 2021, see {doc}`/references`), from its
Zenodo record ([doi:10.5281/zenodo.4708099](https://doi.org/10.5281/zenodo.4708099); MIT
licence, Copyright (C) 2021 Robert Gamache, shipped as `absrbrs/LICENSES/TIPS-2021-MIT.txt`).
The record has one `.QTpy` file per isotopologue, with Q at 1 K steps from 1 K to the table's
maximum (1000 to 9000 K). absrbrs has 181 isotopologues, 876 441 values in all.

In Python they are {func}`absrbrs.partition.tips2021`.

`tools/tips2021.py` builds `absrbrs/data/tips2021.bin` from the record offline, and the result
is committed with its manifest, `absrbrs/data/tips2021.json`. Nothing is unpickled at run
time. The converter:

- checks each file's MD5 against the Zenodo record;
- refuses any pickle opcode beyond those of a plain dict of strings, and loads the file with an
  unpickler that refuses all globals;
- stores Q as f32 (see "Accuracy" below);
- cross-checks every stored value against the Fortran tables of `BD_TIPS_2021.FOR` (in
  `TIPS_2021.zip`, same record, MD5-checked), which give Q to 5 significant figures at 1 to
  20 K by 1 K and 22 to 9000 K by 2 K. Each of the 440 011 shared nodes agrees to within half a
  unit in the 5th significant figure plus the f32 rounding; the worst is at 0.9988 of that
  bound.

### Table ranges

Each table covers its own [T_first, T_last] (`Tips2021::temperature_range` in Rust); a
temperature outside it, or NaN, is an `OutsideTable` error (`OutsideTableError` in Python),
and an isotopologue without a table is `MissingPartitionFunction`. Nothing is extrapolated.

Most tables start at 1 K. Three are trimmed at the low end:

- **H₂S 31_2 and 31_3 start at 30 K.** Their values are unphysical below that: Q(1 K) is
  negative, and the ratios Q(31_2)/Q(31_1) and Q(31_3)/Q(31_1), which should be close to the
  ratios of the nuclear-spin degeneracies (about 1 and 4), are 0.44 and 1.77 at 2 K, 1.43 and
  5.72 at 6 K and 1.03 and 4.11 at 20 K. From 30 K to about 300 K they are within 0.6% of
  them (at most 0.58% and 0.46%). Above that they drift further, to 1.5% at 1000 K and 5.5% at
  4000 K, as expected physically: the isotopologues' rovibrational levels differ.
- **CH₃ 57_1 starts at 2 K,** because its Q(1 K) is 0. Its low-temperature values carry few
  digits (Q(2 K) = 0.000416, about 10⁻³ relative rounding).

`TIPS_2021.py` gives 3010 K as the maximum for CF₄ 42_1, but its table runs, smoothly, to
5000 K, and absrbrs uses all of it.

### Interpolation

Between the nodes, absrbrs interpolates ln Q in T with cubic Hermite polynomials. The node
slopes are the second-order differences of numpy's `gradient` with `edge_order=2`: at an
interior node, the derivative of the quadratic through it and its two neighbours; at an end
node, the derivative of the quadratic through the end three nodes, taken at the end. (On the
uniform 1 K grid these are the central differences and the one-sided second-order ends.)
d ln Q/dT is the derivative of the interpolant itself, so Q and d ln Q/dT are exactly
consistent, and Q = exp(ln Q) is positive everywhere. At the nodes the values are the stored
values.

Interpolating ln Q rather than Q suits its shape, since Q spans orders of magnitude over a
table, and gives d ln Q/dT, the quantity the Jacobians need, directly.

### Accuracy

The interpolation figures below compare absrbrs at the half-integer temperatures of every
table with two other interpolants of the same 1 K nodes: linear interpolation in Q, which is
what the record's own `TIPS_2021.py` does, and 6-point Lagrange interpolation of ln Q, a
higher-order interpolant whose distance from absrbrs estimates the uncertainty between the
nodes. They are reproduced by
`uv run --locked --directory validation python -m absrbrs_validation.tips_interpolation`.

- Q is stored as f32, which rounds by at most 6 × 10⁻⁸ relative, comparable to the 6 to 9
  significant figures of the source.
- Between the nodes, from 50 K up and more than 4 K from the kinks below, absrbrs agrees with
  6-point Lagrange interpolation of ln Q to a median of 3 × 10⁻⁸ and at most 4 × 10⁻⁷ in ln Q;
  d ln Q/dT carries about 10⁻⁵ relative noise from the data.
- Below about 20 K, where ln Q curves sharply, values between the nodes are uncertain by
  10⁻² to 2 × 10⁻¹: absrbrs differs from 6-point Lagrange interpolation of ln Q by up to
  10⁻² (NO⁺ 36_1), and from linear interpolation in Q by up to 0.16 (0.58 for CH₃ 57_1).
- Kinks in the source tables (SO₃ 47_1 at 650 K, HCN 23_1 and 23_2 at 1000 K) distort
  d ln Q/dT within ±1 K of them; within 4 K of them absrbrs differs from 6-point Lagrange
  interpolation of ln Q by up to 5 × 10⁻⁵.
- At non-integer temperatures absrbrs differs from `TIPS_2021.py`'s linear interpolation in Q
  by a median of 4 × 10⁻⁶ and at most 2.7 × 10⁻⁵ (SF₆ 30_1) at 150 to 350 K, by at most
  6.4 × 10⁻⁵ from 50 K (35_1), and by 1.4 × 10⁻⁴ at the SO₃ kink.

## User tables

`TabulatedPartition` (Rust only so far) holds tables of Q(T) supplied by the user, one per
isotopologue; inserting a table replaces any earlier one for the isotopologue. A table needs
at least two points, with T finite, positive and strictly increasing (and ln T strictly
increasing in floating point), and Q finite and positive; the grid need not be uniform.
Outside a table's range evaluation is again `OutsideTable`.

User tables are interpolated in ln Q against **ln T**, not T, with the same cubic Hermite
interpolant and the same node slopes (on the nonuniform ln T grid); d ln Q/dT is
(d ln Q/d ln T)/T. That is exact for Q ∝ Tⁿ, a good model for coarse tables; for a two-point
table the slopes are the secant, which is exact for a power law through both points.
TIPS-2021, with its 1 K grid, interpolates in T instead.

## Why HAPI's tables are not the reference

HAPI 1.3.0.0 (Kochanov et al. 2016) embeds TIPS-2021 at 10 K steps (1, 10, 20, … K), and
`partitionSum(..., version=2021)` (through `BD_TIPS_2021_PYTHON` and `AtoB`) interpolates
them with 3- and 4-point Lagrange polynomials in Q. At integer temperatures absrbrs returns the
1 K source, so comparing the two at the integer temperatures between HAPI's nodes measures
HAPI's error against the source (same script as above):

- up to 1.1 × 10⁻⁵ at 200 to 300 K (SF₆ 30_1, the worst of the 144 isotopologues in both);
- under 10⁻⁶ at 150 to 350 K for the main atmospheric species (molecules 1 to 7,
  isotopologue 1; the worst is O₃ at 9.5 × 10⁻⁷);
- up to about 100% below 20 K (NO⁺ 36_1).

So absrbrs takes the 1 K source as the reference, not HAPI's `partitionSum`.

`validation/tests/test_tips.py` checks both cases for the 144 isotopologues that both have:

- At HAPI's nodes, from 10 K to 1000 K (skipping 10 and 20 K for 31_2 and 31_3), the two agree
  to the rounding of HAPI's tables, which hold Q to 7 significant figures: each node is within
  half a unit in the 7th significant figure of HAPI's value plus the f32 rounding of absrbrs's.
  The largest relative difference is 5.5 × 10⁻⁷, from HAPI's rounding.
- Between HAPI's nodes, at every non-node integer temperature from 150 K to the end of HAPI's
  table, HAPI is within 1.2 × 10⁻⁵ of the source (worst 1.1 × 10⁻⁵, SF₆ 30_1 at 265 K), except
  within 20 K of the three kinks, where its stencil spans the slope discontinuity and the error
  reaches 2.0 × 10⁻³ (SO₃ 47_1 at 651 K).

HAPI 1.3.0.0 has no TIPS-2021 table for 56_1 or 57_1. Those two, and the 35 tables whose
isotopologues are not in `hapi.ISO`, 37 of the 181 in all, are checked only against the
Fortran tables.
