"""
Supertype for line-search specifications. Extend [`configure_linesearch`](@ref)
and [`perform_linesearch!`](@ref) to implement a method.
"""
abstract type AbstractLineSearchMethod end

"""
    LineSearchResult(success, alpha, reason, message, nfev, ndev,
                     iterations, x_trial, F_trial)

Search result with these fields:

- `success`, `alpha`: acceptance flag and step; on failure, the step is diagnostic only.
- `reason`, `message`: termination symbol and optional error detail.
- `nfev`, `ndev`: value and directional-derivative requests.
- `iterations`: step reductions for BackTracking.
- `x_trial`, `F_trial`: borrowed work buffers; accepted data only on success.

Counts reset per search and increment before callbacks, including reported
failures. A joint callback increments both counters once, even if one output
is discarded. Reused data do not count again; work/cache inside callbacks is
not tracked. Their sum is not a count of callback calls or computational cost.
Copy accepted arrays before reuse; see [Using a line search](@ref).
"""
struct LineSearchResult{T<:AbstractFloat,X<:AbstractVector,F<:AbstractVector}
    success::Bool
    alpha::T
    reason::Symbol
    message::Union{Nothing,String}
    nfev::Int
    ndev::Int
    iterations::Int
    x_trial::X
    F_trial::F
end
