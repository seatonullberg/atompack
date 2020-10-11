# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

* `spacegroup` module.

### Changed

* `atom.Atom` and `bond.Bond` now inherit from `MutableMapping`.
* `crystal.Crystal` and `crystal.UnitCell` are now initialized by spacegroup.


## [0.3.0] - 2020-09-08

### Added

* `bond` module.
* `molecule` module.
* `topology` module.
* `crystal.UnitCell` class with many subclasses.
* Support for `pipenv` installation.
* Additional benchmark tests.

### Changed

* `crystal.Crystal` now inherits from `topology.Topology`.

### Removed

* C extensions
* `structure` module.


## [0.2.1] - 2020-08-02

### Added

* `MANIFEST.in`.

### Changed

* Fix `setup.py` to render long description content properly on PyPI.


## [0.2.0] - 2020-07-11

### Changed

* Replace `AtomCollection` with `Structure`.
* Replace pure Python logic with C extensions for a performance boost.
* Separate unit tests and benchmarks into `tests` and `benches` respectively.

### Removed

* `BaseCrystallographyError` and all its subclasses.


## [0.1.0] - 2020-06-03

Initial release.