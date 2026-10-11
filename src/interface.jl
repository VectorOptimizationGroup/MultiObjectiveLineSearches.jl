"""
    configure_linesearch(spec, ::Type{T}) where T<:AbstractFloat

Prepare a line search for vectors with element type `T`. Reuse the returned
configuration across iterations. Invalid options raise `ArgumentError`.
"""
function configure_linesearch end

"""
    perform_linesearch!(search, evaluate!, x, d, Fx, slopes,
                        x_trial, F_trial; factor)

Run a configured line search and return a [`LineSearchResult`](@ref).
On success, the trial buffers hold the accepted point and objectives;
on failure, the caller retains its previous solution. Inputs are preserved.
See [Using a line search](@ref) for evaluator and buffer requirements.
"""
function perform_linesearch! end
