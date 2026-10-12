# MultiObjectiveLineSearches.jl

MultiObjectiveLineSearches.jl provides line searches for multiobjective
optimization. Given a point, a search direction and the objective values at
the point, a line search looks for a step size along the direction that
satisfies an acceptance condition on all objectives. The available method,
`BackTracking`, uses an Armijo-type sufficient decrease condition.
The package has no runtime dependencies and requires Julia 1.11 or later.

## Installation

The package is not registered. From a separate Julia environment, add a local
checkout:

```julia
import Pkg
Pkg.develop(path = "/path/to/MultiObjectiveLineSearches.jl")
```

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

The documentation contains a usage guide, the mathematical formulation of each
line search with its literature sources, the API reference, and the
bibliography. Instructions for building it locally are in
[`docs/README.md`](docs/README.md).

Run the tests from the package directory:

```sh
julia --startup-file=no --project=. test/runtests.jl
```

## License

MultiObjectiveLineSearches.jl is distributed under the [MIT License](LICENSE).
