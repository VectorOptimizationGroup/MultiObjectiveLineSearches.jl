"""
Multiobjective line searches with componentwise sufficient decrease.

Use `configure_linesearch` to prepare a search and `perform_linesearch!`
to find a step along a supplied direction.
"""
module MultiObjectiveLineSearches

export AbstractLineSearchMethod, InitialStatic
export LineSearchResult, configure_linesearch, perform_linesearch!

include("types.jl")
include("initial_step.jl")
include("interface.jl")

end
