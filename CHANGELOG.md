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
