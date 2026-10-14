# API reference

This is the Python API. The Rust API reference will be on [docs.rs](https://docs.rs/absrbrs)
once the crate is published. Both follow the {doc}`conventions`.

## `absrbrs`

```{eval-rst}
.. automodule:: absrbrs
   :no-members:

.. py:data:: __version__
   :type: str

   The installed package's version: the version of the Rust workspace, in PEP 440 form (the
   pre-release ``0.1.0-alpha.1``, for example, is ``0.1.0a1``).
```

The package also re-exports {class}`absrbrs.LineList` and the exceptions of {mod}`absrbrs.errors`,
and provides {mod}`absrbrs.partition` and {mod}`absrbrs.special`.

## `absrbrs.LineList`

```{eval-rst}
.. automodule:: absrbrs.lines
   :no-members:

.. autoclass:: absrbrs.LineList
   :members:
```

## `absrbrs.partition`

```{eval-rst}
.. automodule:: absrbrs.partition
   :members:
```

## `absrbrs.special`

```{eval-rst}
.. automodule:: absrbrs.special
   :members:
```

## `absrbrs.errors`

```{eval-rst}
.. automodule:: absrbrs.errors
   :members:
   :show-inheritance:
```
