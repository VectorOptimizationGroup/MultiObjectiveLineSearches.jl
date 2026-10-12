"""Supertype for multiobjective line-search methods."""
abstract type AbstractLineSearchMethod end

"""
    LineSearchResult(success, alpha, reason, message, nfev, iterations)

Result of a line search. `alpha` is the accepted step on success; on failure,
it is the last candidate step, which may not have been evaluated. `nfev` counts
evaluator calls, including one that reports an error, and `iterations` counts
step reductions. On success, the accepted point and objective values are in the
trial buffers passed to [`perform_linesearch!`](@ref).

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
