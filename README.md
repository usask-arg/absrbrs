# absrbrs

Molecular absorption cross sections, written in Rust with a Python interface. The name reads as
both "absorbers" and "absrb-rs".

Every cross section comes with its Jacobian: with respect to temperature, pressure and volume
mixing ratio, and on request with respect to the spectroscopic line parameters. The line-shape
kernels have bounded, validated accuracy. absrbrs will consolidate the absorption code of
[sasktran2](https://github.com/usask-arg/sasktran2), which will depend on it.

**Status:** pre-alpha. Phase 0 (foundations) is in progress; there is no usable cross-section
API yet. See the
[design specification](https://github.com/usask-arg/absrbrs/blob/main/docs/specs/2026-10-14-absrbrs-design.md)
and [CONTRIBUTING.md](https://github.com/usask-arg/absrbrs/blob/main/CONTRIBUTING.md).

License: MIT.
