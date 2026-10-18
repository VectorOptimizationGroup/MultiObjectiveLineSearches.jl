"""
Multiobjective line searches with reusable preparation and in-place execution.
See [`configure_linesearch`](@ref), [`perform_linesearch!`](@ref) and
[`BackTracking`](@ref).
"""
module MultiObjectiveLineSearches

export AbstractLineSearchMethod, BackTracking, InitialStatic
export LineSearchResult, configure_linesearch, perform_linesearch!

include("types.jl")
include("initial_step.jl")
include("interface.jl")
include("methods/backtracking.jl")

end
