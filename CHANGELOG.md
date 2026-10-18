# Changelog

Notable user-visible changes to MultiObjectiveLineSearches.jl are documented
here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and releases follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- **Breaking:** prepare with `configure_linesearch(spec, evaluate!, x, Fx)`;
  execute with `perform_linesearch!(search, x, d, Fx; slopes, factor, maxstep)`.
  Evaluators and work buffers now belong to the prepared method; old signatures
  have no compatibility overloads.
- Results borrow trial buffers and add the derivative counter `ndev`; a joint
  value/derivative callback increments both `nfev` and `ndev`.
- Per-call `maxstep` restricts the step conservatively; an incompatible bound
  returns `:step_limit` without evaluation.
- Local methods extend the same two functions with their own state and callbacks.

### Added

- `BackTracking`: backtracking for an Armijo-type sufficient decrease
  condition in every objective, with safeguarded quadratic interpolation and
  configurable `c_1`, `tau_1`, `tau_2`, step limits and reduction budget.
- Preparation through `configure_linesearch` and in-place execution through
  `perform_linesearch!`, with method-owned work buffers.
- `LineSearchResult` with borrowed accepted/scratch arrays, step, termination
  reason, message, value/derivative counts and reduction count.
- Documentation relating the condition to the multiobjective literature, with
  executable examples.

[Unreleased]: https://github.com/VectorOptimizationGroup/MultiObjectiveLineSearches.jl/commits/main
