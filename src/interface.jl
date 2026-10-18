"""
    configure_linesearch(spec, evaluate!, x::AbstractVector, Fx::AbstractVector;
                         derivative! = nothing, value_derivative! = nothing)

Prepare reusable method state and work buffers without evaluating objectives.
Fixed options are converted to `T = eltype(x)`; `eltype(Fx)` must equal `T`,
a concrete floating type. Both prototypes must be finite and nonempty.

`evaluate!(F, y)` fills objective values and returns `nothing` or an error
string. Optional derivative callbacks are described in [Extending the protocol](@ref).
See [Using a line search](@ref) for reuse and input requirements.
"""
function configure_linesearch end

"""
    perform_linesearch!(search, x::AbstractVector, d::AbstractVector,
                        Fx::AbstractVector; maxstep=Inf, method_data...)

Run a prepared search and return a [`LineSearchResult`](@ref).
`maxstep` is a nonnegative upper bound on the step (`Inf` adds no restriction);
other keywords depend on the method.

On success, copy the borrowed result arrays into the consumer state before
another call. Inputs must not alias those buffers. See [Using a line search](@ref)
for validation, acceptance and failure handling.
"""
function perform_linesearch! end
