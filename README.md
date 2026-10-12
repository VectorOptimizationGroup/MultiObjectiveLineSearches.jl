# MultiObjectiveLineSearches.jl

[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://vectoroptimizationgroup.github.io/MultiObjectiveLineSearches.jl/dev/)
[![CI](https://github.com/VectorOptimizationGroup/MultiObjectiveLineSearches.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/VectorOptimizationGroup/MultiObjectiveLineSearches.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/VectorOptimizationGroup/MultiObjectiveLineSearches.jl/graph/badge.svg?branch=main)](https://codecov.io/gh/VectorOptimizationGroup/MultiObjectiveLineSearches.jl)

MultiObjectiveLineSearches.jl provides line searches for multiobjective
optimization. Given a point, a search direction and the objective values at
the point, a line search looks for a step size along the direction that
satisfies an acceptance condition on all objectives. The available method,
`BackTracking`, uses an Armijo-type sufficient decrease condition.
The package has no runtime dependencies and requires Julia 1.11 or later.

## Installation

The package is not registered. Add it from its GitHub repository:

```julia
import Pkg
Pkg.add(url = "https://github.com/VectorOptimizationGroup/MultiObjectiveLineSearches.jl")
```

See the [installation guide](https://vectoroptimizationgroup.github.io/MultiObjectiveLineSearches.jl/dev/installation/)
for project environments and local checkouts.

## Quick Example

```julia
using MultiObjectiveLineSearches

function objectives!(values, x)
    values[1] = (x[1] - 1)^2
    values[2] = (x[1] - 1)^2 + 1
    return nothing
end

x = [0.0]; d = [3.0]; Fx = [1.0, 2.0]; slopes = [-6.0, -6.0]

search = configure_linesearch(BackTracking(), Float64)
x_trial = similar(x); F_trial = similar(Fx)

result = perform_linesearch!(search, objectives!, x, d, Fx, slopes,
                             x_trial, F_trial; factor = maximum(slopes))
result.success  # true; x_trial and F_trial hold the accepted point and values
```

## Documentation

The [MultiObjectiveLineSearches.jl documentation](https://vectoroptimizationgroup.github.io/MultiObjectiveLineSearches.jl/dev/)
contains the quick start, the usage guide, the mathematical formulation of
each line search with its literature sources, the API reference, and the
bibliography. Instructions for building it locally are in
[`docs/README.md`](docs/README.md).

Run the tests from the package directory:

```sh
julia --project -e 'using Pkg; Pkg.test()'
```

## Highlights

- Backtracking for an Armijo-type sufficient decrease condition in every
  objective, with safeguarded quadratic interpolation.
- Typed, reusable configurations and in-place searches with caller-owned
  buffers.
- Structured diagnostics: step, termination reason, evaluation and reduction
  counts.
- Formulations related to the multiobjective literature, with executable
  examples.

## Citation

If you use MultiObjectiveLineSearches.jl in research, please cite the software
using the metadata in [`CITATION.cff`](CITATION.cff).

## License

MultiObjectiveLineSearches.jl is distributed under the [MIT License](LICENSE).
