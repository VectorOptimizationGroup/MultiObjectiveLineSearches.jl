"""
Multiobjective line searches with componentwise sufficient decrease.

Use `BackTracking` to configure Armijo backtracking and `perform_linesearch!`
to find a step along a supplied direction.
"""
module MultiObjectiveLineSearches

export AbstractLineSearchMethod, BackTracking, InitialStatic
export LineSearchResult, configure_linesearch, perform_linesearch!

include("types.jl")
include("initial_step.jl")
include("interface.jl")
include("methods/backtracking.jl")

end
