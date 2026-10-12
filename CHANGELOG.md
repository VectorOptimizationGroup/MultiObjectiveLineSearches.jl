# Changelog

Notable user-visible changes to MultiObjectiveLineSearches.jl are documented
here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `BackTracking`: backtracking for an Armijo-type sufficient decrease
  condition in every objective, with safeguarded quadratic interpolation and
  configurable `c_1`, `tau_1`, `tau_2`, step limits and reduction budget.
- Reusable typed configurations through `configure_linesearch`, and in-place
  searches through `perform_linesearch!` with caller-owned buffers.
- `LineSearchResult` with the step, termination reason, message, evaluation
  count and reduction count.
- Documentation relating the condition to the multiobjective literature, with
  executable examples.

[Unreleased]: https://github.com/VectorOptimizationGroup/MultiObjectiveLineSearches.jl/commits/main
