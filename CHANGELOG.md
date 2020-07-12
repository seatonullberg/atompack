# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2020-07-11

### Changed

* Replace `AtomCollection` with `Structure`.
* Replace pure Python logic with C extensions for a performance boost.
* Separate unit tests and benchmarks into `tests` amd `benches` respectively.

### Removed

* `BaseCrystallographyError` and all its subclasses.


## [0.1.0] - 2020-06-03

Initial release.