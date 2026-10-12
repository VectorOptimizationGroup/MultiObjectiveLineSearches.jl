# Installation

MultiObjectiveLineSearches.jl is not registered. The installation commands
below are illustrative: they modify the active environment or use a
placeholder path.

## Install from GitHub

Add MultiObjectiveLineSearches.jl to the active Julia environment directly
from its GitHub repository:

```julia
import Pkg
Pkg.add(url = "https://github.com/VectorOptimizationGroup/MultiObjectiveLineSearches.jl")
```

Using a project-specific environment is recommended so that the package
version is recorded alongside the rest of the project's dependencies. See the
Julia `Pkg` documentation for environment creation and activation. Julia 1.11
or later is required.

## Install from a local checkout

Clone the repository when you need a local copy of the source:

```bash
git clone https://github.com/VectorOptimizationGroup/MultiObjectiveLineSearches.jl.git
```

To install that checkout into the active Julia environment, provide its local
path:

```julia
import Pkg
Pkg.add(path = "/absolute/path/to/MultiObjectiveLineSearches.jl")
```

Use `Pkg.develop(path = ...)` instead when local source edits should be visible
immediately from the active environment.

## Verify the installation

```jldoctest installation_verify
julia> using MultiObjectiveLineSearches

julia> search = configure_linesearch(BackTracking(), Float64);

julia> search.c_1, search.tau_1, search.tau_2
(0.0001, 0.05, 0.95)
```

Continue with the [Quick Start](@ref) for a first line search.
