# Accuracy of the Faddeeva kernels

Every CI run checks both tiers against Arb values computed during the run (`validation/`): about
148 000 grid points, which include both sides of every Fast-tier boundary and of Johnson's branch
thresholds along and near the real axis, and 200 000 random points over the validated domain. The errors use the metrics of
{doc}`/theory/faddeeva`. The Rust tests also check both tiers against the committed Arb grid,
and the Fast tier against the Accurate tier on 20 000 random points and around every region
boundary.

The latest measurement, written by `uv run --directory validation python -m absrbrs_validation
faddeeva`:

```{include} faddeeva-accuracy-table.md
```

The full report, with where each maximum occurs, is
{download}`faddeeva-accuracy.json <data/faddeeva-accuracy.json>`.
