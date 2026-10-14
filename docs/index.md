# absrbrs

Molecular absorption cross sections with first-class linearization, written in Rust with a
Python interface.

```{note}
absrbrs is in early development (phase 0, foundations) and has no cross-section API yet. It is
not published on PyPI or crates.io; see {doc}`installation` to build it.
```

absrbrs computes line-by-line absorption cross sections from spectroscopic line lists. Every
cross section comes with its Jacobian: with respect to temperature, pressure and volume mixing
ratio, and on request with respect to the line parameters. The line-shape kernels have
bounded accuracy, validated against arbitrary-precision references.

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
conventions
```

```{toctree}
:maxdepth: 2
:caption: Validation

validation/faddeeva-prototype
```

```{toctree}
:maxdepth: 2
:caption: Reference

api
references
```
