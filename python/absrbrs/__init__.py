"""absrbrs: molecular absorption cross sections with first-class linearization."""

from importlib.metadata import PackageNotFoundError as _PackageNotFoundError
from importlib.metadata import version as _version

from absrbrs import _absrbrs
from absrbrs.errors import (
    AbsrbrsError,
    InvalidInputError,
    MissingPartitionFunctionError,
    OutsideTableError,
    ParseError,
)

__all__ = [
    "AbsrbrsError",
    "InvalidInputError",
    "MissingPartitionFunctionError",
    "OutsideTableError",
    "ParseError",
    "__version__",
]

# The distribution's version, in PEP 440 form: maturin derives it from the Cargo workspace
# version, which the extension reports as `absrbrs._absrbrs.__version__` in SemVer form (a
# pre-release 0.1.0-alpha.1 there is 0.1.0a1 here). Without distribution metadata (the package
# imported from a source tree with a built extension), the extension's version is used.
try:
    __version__: str = _version("absrbrs")
except _PackageNotFoundError:
    __version__ = _absrbrs.__version__
