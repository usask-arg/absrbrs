import pickle
import sys

import pytest

import absrbrs
from absrbrs._absrbrs import _testing


def test_hierarchy():
    assert issubclass(absrbrs.AbsrbrsError, Exception)
    for cls, std in (
        (absrbrs.InvalidInputError, ValueError),
        (absrbrs.MissingPartitionFunctionError, LookupError),
        (absrbrs.OutsideTableError, ValueError),
        (absrbrs.ParseError, ValueError),
    ):
        assert issubclass(cls, absrbrs.AbsrbrsError)
        assert issubclass(cls, std)
    assert absrbrs.errors.AbsrbrsError is absrbrs.AbsrbrsError


@pytest.mark.parametrize(
    ("kwargs", "cls", "message", "attributes"),
    [
        (
            {"variant": "InvalidInput", "message": "bad"},
            absrbrs.InvalidInputError,
            "invalid input: bad",
            {},
        ),
        (
            {"variant": "MissingPartitionFunction", "molecule": 2, "isotopologue": 11},
            absrbrs.MissingPartitionFunctionError,
            "no partition function for molecule 2, isotopologue 11",
            {"molecule": 2, "isotopologue": 11},
        ),
        (
            {"variant": "OutsideTable", "molecule": 2, "isotopologue": 1, "temperature": 6000.5},
            absrbrs.OutsideTableError,
            "temperature 6000.5 K is outside the partition-function table of molecule 2, "
            "isotopologue 1",
            {"molecule": 2, "isotopologue": 1, "temperature": 6000.5},
        ),
        (
            {"variant": "Parse", "path": "CO2.par", "line_number": 17, "message": "bad width"},
            absrbrs.ParseError,
            "CO2.par, line 17: bad width",
            {"path": "CO2.par", "line_number": 17},
        ),
    ],
)
def test_core_errors_map_to_their_classes(kwargs, cls, message, attributes):
    with pytest.raises(cls) as info:
        _testing.raise_core_error(**kwargs)
    assert type(info.value) is cls
    assert str(info.value) == message
    for name, value in attributes.items():
        assert getattr(info.value, name) == value


@pytest.mark.parametrize(
    ("kind", "cls"),
    [("NotFound", FileNotFoundError), ("PermissionDenied", PermissionError), ("Other", OSError)],
)
def test_io_errors_are_oserrors_of_their_kind(kind, cls):
    with pytest.raises(OSError) as info:
        _testing.raise_core_error("Io", path="CO2.par", kind=kind, message="failure")
    assert type(info.value) is cls
    assert str(info.value) == "cannot read CO2.par: failure"


def test_unknown_variant_and_kind_are_rejected():
    with pytest.raises(ValueError, match="unknown absrbrs::Error variant"):
        _testing.raise_core_error("Nope")
    with pytest.raises(ValueError, match="unknown std::io::ErrorKind"):
        _testing.raise_core_error("Io", kind="Nope")


def test_falls_back_to_runtime_error_when_the_classes_cannot_be_imported(monkeypatch):
    monkeypatch.setitem(sys.modules, "absrbrs.errors", None)
    with pytest.raises(RuntimeError, match="invalid input: bad") as info:
        _testing.raise_core_error("InvalidInput", message="bad")
    assert isinstance(info.value.__cause__, ImportError)


@pytest.mark.parametrize(
    "err",
    [
        absrbrs.MissingPartitionFunctionError("m", 2, 11),
        absrbrs.OutsideTableError("m", 2, 1, 7000.0),
        absrbrs.ParseError("m", "CO2.par", 3),
    ],
)
def test_errors_pickle_with_their_attributes(err):
    back = pickle.loads(pickle.dumps(err))
    assert type(back) is type(err)
    assert str(back) == str(err)
    assert vars(back) == vars(err)
