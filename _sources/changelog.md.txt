# Changelog

## Version 0.4.21

- Optimize `intersect` operation on large number of genomic regions
- Add a `fast_combine_granges` method that only concatenates seqnames and intervals.
- Use `np.asarray` as much as possible for conversion
- Update tests and documentation

## Version 0.4.20

Updated dependencies (especially `IRanges`), to be compatible with NumPy's 2.0 release.

## Version 0.4.19

- Fix an issue related to combining `GenomicRanges` or `GenomicRangesList` objects containing different metadata columns. This switches the operation from `combine_rows` to `relaxed_combine_rows`.
- Set NumPy to <2.0.0 before the migration to 2.0 release.
- Updated tests

## Version 0.4.17

Initialize `GenomicRanges` from a polars `DataFrame`. Also coerce a GRanges object
to polars.

## Version 0.4.9 - 0.4.16

Bring `SeqInfo` upto speed with the rest of the class implementations. Added subset and pretty print functionality.

Fix minor bugs and update documentation.

## Version 0.4.0 to 0.4.8

This is a complete rewrite of both these classes following the functional paradigm from our [developer notes](https://github.com/BiocPy/developer_guide#use-functional-discipline).

`GenomicRanges` is now closer to Bioconductor's GenomicRanges class both in the design and implementation.

The package does not rely on pandas anymore. While we try to provide backwards compatibility to construct a GenomicRanges object from a pandas dataframe using the `from_pandas` method, please note that the default constructor to genomic ranges does not accept a pandas data frame anymore!

Most range based methods have been reimplemented and the heavy lifting is done in the [IRanges package](https://github.com/BiocPy/IRanges) for interval operations. The package indirectly depends on [NCLS](https://github.com/pyranges/ncls) interval tree data structure to perform search and overlap operations.

Tests, documentation and readme has been updated to reflect these changes.

## Version 0.3.0

This release migrates the package to a more palatable Google's Python style guide. A major modification to the package is with casing, all `camelCase` methods, functions and parameters are now `snake_case`.

In addition, docstrings and documentation has been updated to use sphinx's features of linking objects to their types. Sphinx now also documents private and special dunder methods (e.g. `__getitem__`, `__copy__` etc). Intersphinx has been updated to link to references from dependent packages.

Configuration for flake8, ruff and black has been added to pyproject.toml and setup.cfg to be less annoying.

Finally, pyscaffold has been updated to use "myst-parser" as the markdown compiler instead of recommonmark. As part of the pyscaffold setup, one may use pre-commits to run some of the routine tasks of linting and formatting before every commit. While this is sometimes annoying and can be ignored with `--no-verify`, it brings some consistency to the code base.

## Version 0.2

- Now uses BiocFrame as the underlying class
- Implement interval based operations
- update documentation, readme
- comprehensive tests

## Version 0.1

- Create base class for GenomicRanges
- Import into GenomicRanges from pandas/GTF/UCSC
- tests
- documentation
