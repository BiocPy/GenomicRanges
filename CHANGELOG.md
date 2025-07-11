# Changelog

## Version 0.7.0 - 0.7.1

- Changes to switch to LTLA/nclist-cpp in the iranges package for overlap and search operations.
- Improve performance of search operations, bump version of iranges to 0.5.2.
- Optimize group by operations that `GenomicRanges` uses internally for the inter-range operations.

## Version 0.6.2 - 0.6.3

- Implement biocutil's `extract_row_names` generic.
- Fix when SeqInfo's attributes contain lists with numpy elements.

## Version 0.6.0 - 0.6.1

An rewrite of the package to use the new and improve IRanges packages (>= 0.4.2)

- More consistent results across all methods compared to R/Bioconductor implementations.
- Similar to IRanges, search and overlap operations may return a hits like object.
- Nearest returns slightly different matches but since the select="arbitrary", its probably ok.
- More robust testing of the combination of inputs and parameter choices.
- Updated the tests, docstrings

## Version 0.5.2

- Restrict IRanges to the last compatible version before the migration.

## Version 0.5.1

- Fixed an issue with numpy arrays as slice arguments. Code now uses Biocutils's subset functions to perform these operations.
- Rename GitHub actions for consistency with the rest of the packages.

## Version 0.5.0

- chore: Remove Python 3.8 (EOL)
- precommit: Replace docformatter with ruff's formatter

## Version 0.4.32 - 0.4.33

- Bump IRanges package version to fix coercion issues to pandas.
- Remove reverse mapping in iranges in `reduce` operation.
- Fixes issue with combine merging sequence names without properly using the accessor methods.

## Version 0.4.27 - 0.4.31

- Implement `subtract` method, add tests.
- Use accessor methods to access properties especially `get_seqnames()`
- Modify search and overlap methods for strand-awareness.
- Choose appropriate NumPy dtype for sequences.
- Update tests and documentation.

## Version 0.4.26

- Expose the `generic_accessor` method used internally by `GenomicRangesList` to call functions from the underlying GenomicRanges for each element.
- Add class method to initialize `GenomicRangesList` from dictionary.
- Update tests

## Version 0.4.25

- Method to split `GenomicRanges` by a list of groups.
- Coerce `GenomicRangesList` to `GenomicRanges`.
- Add tests and documentation.


## Version 0.4.21 - 0.4.24

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

Most range based methods have been reimplemented and the heavy lifting is done in the [IRanges package](https://github.com/BiocPy/IRanges) for interval operations. The package depends on [nclist-cpp](https://github.com/LTLA/nclist-cpp) data structure to perform search and overlap operations.

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
