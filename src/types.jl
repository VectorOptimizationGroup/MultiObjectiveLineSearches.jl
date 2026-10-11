"""Supertype for multiobjective line-search methods."""
abstract type AbstractLineSearchMethod end

"""
    LineSearchResult(success, alpha, reason, message, nfev, iterations)

Result of a line search. `alpha` is the accepted step on success, or the last
candidate on failure; a rejected candidate may not have been evaluated.
`nfev` counts evaluator calls, including reported errors, and `iterations`
counts step reductions. Accepted vectors are returned through the trial buffers.

`reason` is `:success`, `:max_iterations`, `:minimum_step`,
`:evaluation_failure`, or `:indistinguishable_point`. `message` contains
additional detail, or `nothing`.
"""
struct LineSearchResult{T<:AbstractFloat}
    success::Bool
    alpha::T
    reason::Symbol
    message::Union{Nothing,String}
    nfev::Int
    iterations::Int
end
