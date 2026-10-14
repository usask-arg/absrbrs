"""Registry of reference generators. Each returns ``(precision_bits, records)``."""

from collections.abc import Callable

from absrbrs_reference.generators import faddeeva, weideman

Generator = Callable[[], tuple[int, list[dict]]]

GENERATORS: dict[str, Generator] = {
    "faddeeva": faddeeva.generate,
    "weideman": weideman.generate,
}
