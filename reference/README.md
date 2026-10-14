# absrbrs reference generator

Independent ground truth for absrbrs tests, computed with Arb ball arithmetic through
[python-flint](https://python-flint.readthedocs.io). It does not import absrbrs.

    uv run --directory reference python -m absrbrs_reference generate            # all files
    uv run --directory reference python -m absrbrs_reference generate --only faddeeva
    uv run --directory reference python -m pytest

The record format is that of the sibling project scatrs. Each value is stored as three doubles
`hi`/`lo`/`lo2` (`hi` is the double nearest to the value, `hi + lo` a normalised double-double
and `lo2` the double nearest to the rest, about 160 bits in all), a decimal string `dec` (up to
25 significant digits, fewer when the ball is wider), and a rigorous error bound `err` on
|true − (hi + lo + lo2)|. `label` is `certified` when every source of error is rigorously
bounded, otherwise `high-precision`.

`data/manifest.json` records, for each file, the generator, the working precision, the
python-flint and FLINT versions and the SHA-256 of the file. python-flint is pinned exactly:
CI regenerates every file and fails on any difference from the committed bytes.

| File | Contents |
|---|---|
| `faddeeva.jsonl` | w(z) = exp(−z²) erfc(−iz) and w′(z) on the grid of spec §6.2: 124 values of x from 0 to 1e8 and 74 values of y from 0 to 1e6 |

Generators draw their inputs without the platform's maths library (use `nearest_pow10` from
`absrbrs_reference/records.py` instead of `10.0 ** u`), so the data are the same bytes on every
platform; `tests/test_portable_inputs.py` enforces this.
