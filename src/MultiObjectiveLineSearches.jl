"""
Line searches for multiobjective optimization.

Prepare a method such as `BackTracking` with `configure_linesearch`, then call
`perform_linesearch!` to search for a step size along a given direction.
"""
module MultiObjectiveLineSearches

export AbstractLineSearchMethod, BackTracking, InitialStatic
export LineSearchResult, configure_linesearch, perform_linesearch!

include("types.jl")
include("initial_step.jl")
include("interface.jl")
include("methods/backtracking.jl")

end
