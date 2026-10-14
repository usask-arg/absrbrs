# Changelog

All notable changes to absrbrs are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- Cargo workspace with the `absrbrs` core crate and its error type.
- The `absrbrs` Python package (extension module `absrbrs._absrbrs`) with its exception
  hierarchy.
- Arb reference values of the Faddeeva function and its derivative.
- The phase 0 Faddeeva prototype and its report, which chose the Fast-tier algorithm.
- Third-party licence register (`THIRD_PARTY.md`) with the outcome of each phase 0 licence check.
- Contributor tooling: pre-commit hooks, dependabot configuration and `CONTRIBUTING.md`.
- Documentation site, CI, release and scheduled workflows.
- `absrbrs::faddeeva`: w(z) and its derivative, and the Voigt kernel triplet, at the Fast and
  Accurate tiers (the Accurate tier ports the MIT-licensed Faddeeva Package).
- `absrbrs.special.faddeeva` and `absrbrs.special.voigt_kernel`.
- Validation of both tiers against Arb in CI; benchmarks against sasktran2's J-pole
  approximation.
- `absrbrs::lines`: `LineList`, a structure-of-arrays line list sorted by ν₀, with optional
  fields (n_self, δ_self, δ′_air) and their defaults, `select` and `with_column`; the HITRAN
  isotopologue table (masses and abundances, generated from HAPI's `hapi.ISO`).
- `absrbrs::io`: readers for HITRAN `.par` files and HAPI `.data`/`.header` caches.
- `absrbrs::partition`: the `PartitionFunction` trait; the built-in TIPS-2021 tables
  (181 isotopologues at 1 K steps, MIT-licensed, embedded as f32) and user tables
  (`TabulatedPartition`), both giving Q and a consistent d ln Q/dT by cubic Hermite
  interpolation of ln Q.
- The crate-private line-state layer: the line-shape quantities (A, ν₀′, σ_D, y) of a line at
  a condition, with closed-form derivatives with respect to T, p, VMR and every line field,
  checked against dual numbers and Arb.
- `absrbrs.LineList` (with `from_hitran`, `from_hapi` and `from_arrays`) and
  `absrbrs.partition.tips2021`.
- Validation of TIPS-2021 against HAPI's embedded tables; theory pages on partition functions
  and the line state.
