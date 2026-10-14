import importlib
import importlib.metadata
import tomllib
from importlib.metadata import version
from pathlib import Path

from packaging.version import Version

import absrbrs
from absrbrs import _absrbrs

ROOT = Path(__file__).parents[1]


def test_version_matches_distribution_metadata():
    assert absrbrs.__version__ == version("absrbrs")


def test_version_matches_cargo_workspace():
    # The single version source: [workspace.package] version in the root Cargo.toml. It is
    # SemVer; maturin writes it in normalised PEP 440 form (0.1.0-alpha.1 becomes 0.1.0a1).
    cargo = tomllib.loads((ROOT / "Cargo.toml").read_text("utf-8"))
    assert absrbrs.__version__ == str(Version(cargo["workspace"]["package"]["version"]))


def test_extension_is_built_from_the_installed_version():
    assert Version(_absrbrs.__version__) == Version(absrbrs.__version__)


def test_version_falls_back_to_the_extension_without_metadata(monkeypatch):
    def not_installed(name):
        raise importlib.metadata.PackageNotFoundError(name)

    monkeypatch.setattr(importlib.metadata, "version", not_installed)
    try:
        assert importlib.reload(absrbrs).__version__ == _absrbrs.__version__
    finally:
        monkeypatch.undo()
        importlib.reload(absrbrs)
    assert absrbrs.__version__ == version("absrbrs")


def test_crate_license_is_the_repository_license():
    # absrbrs/LICENSE is a copy, so that the published crate ships the license text.
    assert (ROOT / "absrbrs" / "LICENSE").read_text("utf-8") == (ROOT / "LICENSE").read_text(
        "utf-8"
    )


def test_public_names():
    for name in absrbrs.__all__:
        assert hasattr(absrbrs, name)


def test_test_hooks_live_outside_the_extension_namespace():
    assert set(_absrbrs.__all__) == {"__version__"}
    testing = importlib.import_module("absrbrs._absrbrs._testing")
    assert testing is _absrbrs._testing
    assert set(testing.__all__) == {"raise_core_error"}
