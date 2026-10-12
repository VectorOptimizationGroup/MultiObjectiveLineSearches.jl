# Installation

## Install from a local checkout

MultiObjectiveLineSearches.jl is not registered. To install a local checkout
into the active Julia environment, provide its path. The installation command
below is illustrative: it modifies the active environment and uses a
placeholder path.

```julia
import Pkg
Pkg.add(path = "/absolute/path/to/MultiObjectiveLineSearches.jl")
```

Use `Pkg.develop(path = ...)` instead when local source edits should be visible
immediately from the active environment. Using a project-specific environment
is recommended so that the package version is recorded alongside the rest of
the project's dependencies. Julia 1.11 or later is required.

## Verify the installation

```jldoctest installation_verify
julia> using MultiObjectiveLineSearches

julia> search = configure_linesearch(BackTracking(), Float64);

julia> search.c_1, search.tau_1, search.tau_2
(0.0001, 0.05, 0.95)
```

Continue with the [Quick Start](@ref) for a first line search.
