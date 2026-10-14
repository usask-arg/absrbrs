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
cargo test --workspace
cargo +1.85.0 check --workspace --lib   # the MSRV: rustup toolchain install 1.85.0
cargo package -p absrbrs                # add --allow-dirty while changes are uncommitted
uv run ruff check .
uv run ruff format --check .
uv run pytest
uv run python -m mypy.stubtest absrbrs._absrbrs
uv run --directory reference python -m pytest
uv run --directory reference python -m absrbrs_reference generate   # then: git status reference/data
uv run --directory prototypes/faddeeva python -m pytest
uv run --extra docs sphinx-build -E -W --keep-going -b html docs docs/_build/html
uvx --from actionlint-py==1.7.12.25 actionlint   # after changing .github/workflows
```

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
  - "All checks passed", which succeeds only if every job above did: the one check to
    require.
- `.github/workflows/release.yml`, on pushes to `main`, tags `v*`, and pull requests that touch
  the packaging: wheels for every platform and the sdist, installed and tested on every
  platform with a GitHub runner. A tag also runs the preflight checks and publishes.
- `.github/workflows/docs.yml`: publishes the documentation from `main` to GitHub Pages.
- `.github/workflows/scheduled.yml`, every Monday: the Rust tests on Intel macOS, and a
  regeneration of every reference file on Linux and macOS that fails on any difference from
  the committed data.

## Reference data

Arb reference values live in `reference/data/` and are produced by `reference/` (python-flint,
pinned). Regenerate them with
`uv run --directory reference python -m absrbrs_reference generate [--only NAME]`
and commit the data together with the manifest. Never edit the data by hand. Tolerances in
tests are never loosened without a written justification in the test.

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

One-time repository setup: GitHub Pages with source "GitHub Actions"; a tag ruleset that
restricts `v*` tags to maintainers; Trusted Publishers on PyPI and crates.io for this
repository's `release.yml`; the first crates.io release published by hand, because Trusted
Publishing needs an existing crate.
