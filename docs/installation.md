# Installation

absrbrs is not yet published. To build it from a clone of the repository you need a stable
Rust toolchain (1.85 or later), Python 3.11 or later, and [uv](https://docs.astral.sh/uv/):

```bash
uv sync          # builds the extension and installs absrbrs into .venv
uv run pytest    # the Python tests
```

`pip install .` also works; it builds the extension with maturin.

## Building the documentation

```bash
uv run --extra docs sphinx-build -E -W --keep-going -b html docs docs/_build/html
```

The build downloads the Python and numpy inventories for cross-references. Without network
access, set `ABSRBRS_DOCS_OFFLINE=1`: references to Python and numpy objects are then rendered
as plain text.
