# absrbrs

Molecular absorption cross sections with first-class linearization, written in Rust with a
Python interface.

```{note}
absrbrs is in early development: phase 1a (the Faddeeva kernels and `absrbrs.special`) is complete, and there is still no cross-section API. It is not published on PyPI or crates.io; see {doc}`installation` to build it.
```

absrbrs will compute line-by-line absorption cross sections from spectroscopic line lists.
Every cross section is designed to come with its Jacobian: with respect to temperature,
pressure and volume mixing ratio, and on request with respect to the line parameters. The
line-shape kernels are designed to have bounded accuracy, validated against
arbitrary-precision references.

```{toctree}
:maxdepth: 1
:caption: User guide

installation
conventions
```

```{toctree}
:maxdepth: 1
:caption: Theory

theory/faddeeva
```

```{toctree}
:maxdepth: 1
:caption: Validation

validation/faddeeva-accuracy
validation/faddeeva-prototype
```

```{toctree}
:maxdepth: 1
:caption: Reference

api
references
```
