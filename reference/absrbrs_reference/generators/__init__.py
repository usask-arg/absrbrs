"""Registry of reference generators. Each returns ``(precision_bits, records)``."""

from collections.abc import Callable

from absrbrs_reference.generators import faddeeva, line_state, weideman

Generator = Callable[[], tuple[int, list[dict]]]

GENERATORS: dict[str, Generator] = {
    "faddeeva": faddeeva.generate,
    "line_state": line_state.generate,
    "weideman": weideman.generate,
}
