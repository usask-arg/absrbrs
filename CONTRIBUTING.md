# Contributing to absrbrs

## Setup

```bash
rustup update stable        # CI uses the latest stable Rust (rust-toolchain.toml: "stable")
uv sync                     # Python environment with the dev tools; builds the extension
uv run pre-commit install   # formatting and lint checks on every commit
```

`uv run` rebuilds the extension whenever the Rust sources change, so `uv run pytest` always
tests the current code. Keep stable Rust up to date: CI installs the latest, whose new clippy
lints an older local toolchain does not show.

Run cargo from inside the repository. `.cargo/config.toml` sets
`PYO3_BUILD_EXTENSION_MODULE=1`, which lets cargo build `absrbrs-py` without a suitable Python
interpreter, but cargo reads that file only from the directory it is started in and its
parents.

## Checks (the same ones CI runs)

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo build --workspace --lib
RUSTDOCFLAGS="-D warnings" cargo doc -p absrbrs --no-deps --document-private-items
cargo check --workspace --target s390x-unknown-linux-gnu   # rustup target add s390x-unknown-linux-gnu
cargo test --workspace
cargo test -p absrbrs --bench faddeeva   # runs each benchmark once; timing is not checked
cargo +1.85.0 check --workspace --lib   # the MSRV: rustup toolchain install 1.85.0
cargo package -p absrbrs                # add --allow-dirty while changes are uncommitted
uv run ruff check .
uv run ruff format --check .
uv run pytest
uv run python -m mypy.stubtest absrbrs._absrbrs
uv run --directory reference python -m pytest
uv run --directory reference python -m absrbrs_reference generate   # then: git status reference/data
uv run --directory prototypes/faddeeva python -m pytest
uv run --directory validation python -m pytest   # Faddeeva tiers against Arb; TIPS-2021 against HAPI
uv run --extra docs sphinx-build -E -W --keep-going -b html docs docs/_build/html
uvx --from actionlint-py==1.7.12.25 actionlint   # after changing .github/workflows
```

Wall-clock benchmarks of the Faddeeva kernels: `cargo bench -p absrbrs --bench faddeeva`.

CI also tests the packaged crate on its own, as `cargo test` on the crates.io sources would:
extract `target/package/absrbrs-*.crate` outside the repository and run `cargo test --locked`
there (the reference-backed integration tests are not published). The pre-commit hooks
(whitespace, end-of-file, YAML and TOML checks, Markdown formatting) run locally only, not in
CI.

The documentation build downloads the Python and numpy inventories; see
`docs/installation.md` for building it offline.

## Continuous integration

Every action is pinned to a commit SHA, with its version in a comment. Dependabot updates the
actions, `Cargo.lock` and the `uv.lock` files weekly, in one grouped pull request each, and
only to releases at least a week old.

- `.github/workflows/ci.yml`, on every pull request and push to `main`:
  - lint: rustfmt, clippy, rustdoc, `cargo check` for the big-endian wheel target (s390x;
    32-bit platforms are not supported), the crate package and its tests, ruff and actionlint;
  - the Rust tests on Linux, macOS and Windows, and `cargo check` of the libraries with Rust
    1.85.0 (the MSRV);
  - on each of the three systems, one abi3 wheel built and tested with Python 3.11 and 3.14,
    and on Linux also with Python 3.11 and the oldest supported numpy (1.26);
  - the strict documentation build and the type stubs;
  - the reference generator's tests, a regeneration of every reference file, and the Faddeeva
    prototype's tests;
  - the accuracy of both Faddeeva tiers against Arb values computed at run time, and
    TIPS-2021 against HAPI's embedded tables at HAPI's nodes (`validation/`);
  - "All checks passed", which succeeds only if every job above did: the one check to
    require.
- `.github/workflows/release.yml`, on pushes to `main`, tags `v*`, and pull requests that touch
  the packaging: wheels for every platform and the sdist, installed and tested on every
  platform with a GitHub runner. A tag also checks the tag against the Cargo version and publishes.
- `.github/workflows/docs.yml`: publishes the documentation from `main` to GitHub Pages.
- `.github/workflows/scheduled.yml`, every Monday: the Rust tests on Intel macOS, and a
  regeneration of every reference file on Linux and macOS that fails on any difference from
  the committed data.

## Reference data

After changing dependencies in the root `pyproject.toml`, run `uv lock --directory validation`
and commit `validation/uv.lock`: that lock embeds absrbrs's `requires-dist`, so CI's validation
job fails until it is regenerated. The `validation-lock` pre-commit hook checks it.

Arb reference values live in `reference/data/` and are produced by `reference/` (python-flint,
pinned). Regenerate them with
`uv run --directory reference python -m absrbrs_reference generate [--only NAME]`
and commit the data together with the manifest. Never edit the data by hand. Tolerances in
tests are never loosened without a written justification in the test.

## Local HITRAN line data (opt-in)

HITRAN line data is never committed, shipped or uploaded, including in CI artifacts, so the
tests use synthetic records. One opt-in test reads your own HAPI caches (for example
sasktran2's `database/hitran_lines`), checks every table and prints its line count, timing and
skipped records:

```bash
ABSRBRS_HITRAN_LINES_DIR=/path/to/hitran_lines cargo test -p absrbrs --release --test local_hitran -- --nocapture
```

Without the variable it is skipped; CI never sets it.

## Generated data (`tools/`)

Two uv scripts generate committed data; run them only to regenerate it, and commit the outputs
together:

- `uv run tools/isotopologues.py` writes `absrbrs/src/lines/isotopologue_table.rs` from HAPI's
  `hapi.ISO` (hitran-api, pinned in the script's header). Never edit the generated file.
- `uv run tools/tips2021.py` downloads the TIPS-2021 Zenodo record, checks each file's MD5,
  loads the `.QTpy` pickles with an unpickler that refuses all globals, cross-checks the
  result against the record's Fortran tables, and writes `absrbrs/data/tips2021.bin` and its
  manifest `absrbrs/data/tips2021.json`. `tests/test_data_files.py` checks the two against
  each other.

The interpolation figures quoted for TIPS-2021 (in `docs/theory/partition-functions.md` and the
`Tips2021` rustdoc) are printed by
`uv run --locked --directory validation python -m absrbrs_validation.tips_interpolation`.

## The Faddeeva prototype

`prototypes/faddeeva/` is the phase 0 study that chose the Fast-tier algorithm. Its report is
regenerated with `uv run --directory prototypes/faddeeva python -m faddeeva_prototype report`,
which writes `docs/validation/data/faddeeva_prototype.json` and
`docs/validation/faddeeva-prototype-table.md`. The dense Arb map it uses is cached in
`prototypes/faddeeva/.cache/` (not committed).

## Releasing

1. Set `version` in `[workspace.package]` of `Cargo.toml` (the only place it is set), move the
   `[Unreleased]` entries of `CHANGELOG.md` under the new version, and merge to `main`.
2. Tag the merge commit `v<version>` and push the tag. `release.yml` checks that the tag
   matches the Cargo version and is on `main`, builds and tests every wheel, then publishes to
   PyPI and crates.io with Trusted Publishing (environments `pypi` and `crates-io`, each with
   required reviewers).

One-time repository setup:

- GitHub Pages with source "GitHub Actions", before the first push to `main` (`docs.yml`
  deploys there).
- A tag ruleset that restricts `v*` tags to maintainers.
- Create the GitHub environments `pypi` and `crates-io`, each with required reviewers and a
  deployment rule allowing only `v*` tags, before the first tag: a workflow that references a
  missing environment creates it without protection.
- Trusted Publishers on PyPI and crates.io for this repository's `release.yml`.
- First-release order: Trusted Publishing needs an existing crate, so publish the first crate
  version to crates.io by hand and configure its Trusted Publisher before pushing the first
  tag. That first tag's crates.io job then fails with "already exists", which is expected
  (PyPI still publishes); later tags publish both.
