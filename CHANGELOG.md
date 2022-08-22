# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- R package structure
- roxygen documentation
- unit tests with [testthat](https://testthat.r-lib.org/)
- basic power calculation, variance explained, expected replication
- Poisson binomial test, for assessing replication power distribution

### Changed

- function `calculate_ci` uses basic asymptotic normal theory
- input format is somewhat in flux; will be locked down before release version
- standard error reestimation is entirely disabled, due to the removal of non-linear GLM.
  input standard error is now preserved, taking advantage of the fact that, under the linear
  model, genetic association standard error is independent of regression coefficient

### Removed

- entire existing C++ implementation
- function `ci_log_likelihood` was slowly reimplementing normal theory
- logistic and poisson regression were never exactly functional, so they are removed for now

### Fixed

- MSE estimation for coefficient and CI once again follows actual logic
  - this was not providing performance improvements in the original runs, thus it was hackjob
    overridden; but it should never have been

[//]: # (- Added)
[//]: # (- Changed)
[//]: # (- Deprecated)
[//]: # (- Removed)
[//]: # (- Fixed)
[//]: # (- Security)
