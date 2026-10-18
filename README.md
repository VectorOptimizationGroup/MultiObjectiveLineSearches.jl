# MultiObjectiveLineSearches.jl

[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://vectoroptimizationgroup.github.io/MultiObjectiveLineSearches.jl/dev/)
[![CI](https://github.com/VectorOptimizationGroup/MultiObjectiveLineSearches.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/VectorOptimizationGroup/MultiObjectiveLineSearches.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/VectorOptimizationGroup/MultiObjectiveLineSearches.jl/graph/badge.svg?branch=main)](https://codecov.io/gh/VectorOptimizationGroup/MultiObjectiveLineSearches.jl)

Multiobjective line searches with reusable preparation and in-place execution.
`BackTracking` implements an Armijo-type decrease condition with safeguarded
quadratic interpolation. Requires Julia 1.11 or later; no runtime dependencies.

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

search = configure_linesearch(BackTracking(), objectives!, x, Fx)
result = perform_linesearch!(search, x, d, Fx;
                             slopes, factor=maximum(slopes), maxstep=1)
if result.success
    copyto!(x, result.x_trial)
    copyto!(Fx, result.F_trial)
end
# x == [1.0], Fx == [0.0, 1.0]
```

Prepare once and reuse `search`. Copy accepted results before the next call:
their arrays borrow reusable work buffers.

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

## Citation

If you use MultiObjectiveLineSearches.jl in research, please cite the software
using the metadata in [`CITATION.cff`](CITATION.cff).

## License

MultiObjectiveLineSearches.jl is distributed under the [MIT License](LICENSE).
