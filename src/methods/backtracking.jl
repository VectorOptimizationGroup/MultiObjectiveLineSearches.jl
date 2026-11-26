"""
    BackTracking(; initial_step=InitialStatic(), c_1=1e-4, tau_1=0.05, tau_2=0.95,
                 max_iters=1000, minstep=1e-15, maxstep=1.0, verbose=0)

Armijo-type backtracking with safeguarded quadratic interpolation.
`max_iters` counts reductions; `verbose=1` logs the final diagnostic.
Prepare with [`configure_linesearch`](@ref) and supply `slopes` and `factor`
to each execution; `reference` and `alpha0` are optional per call.

See [BackTracking](@ref backtracking) for options, admissible ranges,
acceptance and reduction rules, and literature sources.
"""
struct BackTracking{I<:InitialStatic,C<:Real,A<:Real,B<:Real,L<:Real,U<:Real} <: AbstractLineSearchMethod
    initial_step::I
    c_1::C
    tau_1::A
    tau_2::B
    max_iters::Int
    minstep::L
    maxstep::U
    verbose::Int
end

function BackTracking(; initial_step::InitialStatic=InitialStatic(),
                      c_1::Real=1e-4, tau_1::Real=0.05, tau_2::Real=0.95,
                      max_iters::Integer=1000, minstep::Real=1e-15,
                      maxstep::Real=1.0, verbose::Integer=0)
    return BackTracking(initial_step, c_1, tau_1, tau_2, Int(max_iters),
                        minstep, maxstep, Int(verbose))
end

struct BacktrackingParams{T<:AbstractFloat}
    initial_step::InitialStatic{T}
    c_1::T
    tau_1::T
    tau_2::T
    max_iters::Int
    minstep::T
    maxstep::T
    verbose::Int
end

struct PreparedBackTracking{T<:AbstractFloat,E,X<:AbstractVector,F<:AbstractVector}
    options::BacktrackingParams{T}
    evaluate!::E
    x_trial::X
    F_trial::F
end

function configure_linesearch(spec::BackTracking, evaluate!::E,
                              x::AbstractVector, Fx::AbstractVector;
                              derivative! = nothing, value_derivative! = nothing) where {E}
    T = eltype(x)
    T <: AbstractFloat && isconcretetype(T) ||
        throw(ArgumentError("execution type must be concrete and floating"))
    eltype(Fx) === T || throw(ArgumentError("point and objectives must have the same element type"))
    !isempty(x) && !isempty(Fx) || throw(ArgumentError("prototypes must be nonempty"))
    all(isfinite, x) && all(isfinite, Fx) || throw(ArgumentError("prototypes must be finite"))
    spec.max_iters >= 0 || throw(ArgumentError("max_iters must be nonnegative"))
    spec.verbose in (0, 1) || throw(ArgumentError("verbose must be 0 or 1"))
    c_1, lo, hi, alpha = T(spec.c_1), T(spec.minstep), T(spec.maxstep), T(spec.initial_step.alpha)
    isfinite(c_1) && zero(T) < c_1 < one(T) ||
        throw(ArgumentError("c_1 must be finite and in (0, 1) after conversion"))
    tau_1, tau_2 = T(spec.tau_1), T(spec.tau_2)
    zero(T) < tau_1 <= T(0.5) <= tau_2 < one(T) ||
        throw(ArgumentError("require 0 < tau_1 <= 1/2 <= tau_2 < 1 after conversion"))
    all(isfinite, (lo, hi, alpha)) && zero(T) < lo < hi && lo <= alpha <= hi ||
        throw(ArgumentError("require finite 0 < minstep < maxstep and minstep <= alpha <= maxstep after conversion"))
    options = BacktrackingParams(InitialStatic(alpha), c_1, tau_1, tau_2, spec.max_iters,
                                 lo, hi, spec.verbose)
    x_trial, F_trial = similar(x), similar(Fx)
    eltype(x_trial) === T && eltype(F_trial) === T ||
        throw(ArgumentError("trial buffers must preserve the execution type"))
    axes(x_trial) == axes(x) && axes(F_trial) == axes(Fx) ||
        throw(DimensionMismatch("trial buffers must preserve prototype axes"))
    for buffer in (x_trial, F_trial), input in (x, Fx)
        Base.mightalias(buffer, input) && throw(ArgumentError("trial buffers must not alias prototypes"))
    end
    Base.mightalias(x_trial, F_trial) && throw(ArgumentError("trial buffers must not alias each other"))
    applicable(evaluate!, F_trial, x_trial) ||
        throw(ArgumentError("evaluate! must accept the prepared value and point buffers"))
    return PreparedBackTracking(options, evaluate!, x_trial, F_trial)
end

"""
    _quadratic_step(alpha, value, slope, value0, tau_1, tau_2)

Return the safeguarded quadratic proposal, or `alpha/2` if it is invalid.
See [Quadratic interpolation and safeguards](@ref quadratic-interpolation).
"""
function _quadratic_step(alpha::T, value::T, slope::T, value0::T, tau_1::T, tau_2::T) where {T}
    if slope < zero(T)
        den = value - value0 - alpha * slope
        if !iszero(den) && isfinite(den)
            candidate = (-slope * alpha^2) / (T(2) * den)
            if isfinite(candidate) && tau_1 * alpha <= candidate <= tau_2 * alpha
                return candidate
            end
        end
    end
    return alpha * T(0.5)
end

function _validate_inputs(::Type{T}, x, d, Fx, slopes, x_trial, F_trial, factor) where {T}
    slopes isa AbstractVector || throw(ArgumentError("slopes must be a vector"))
    inputs = (x, d, Fx, slopes)
    buffers = (x_trial, F_trial)
    all(v -> eltype(v) === T, inputs) ||
        throw(ArgumentError("all vectors must have the configured element type $T"))
    !isempty(x) && !isempty(Fx) || throw(ArgumentError("input vectors must be nonempty"))
    axes(x) == axes(d) == axes(x_trial) ||
        throw(DimensionMismatch("point, direction and trial-point axes must match"))
    axes(Fx) == axes(slopes) == axes(F_trial) ||
        throw(DimensionMismatch("objective, slope and trial-objective axes must match"))
    for buffer in buffers, input in inputs
        Base.mightalias(buffer, input) && throw(ArgumentError("trial buffers must not alias inputs"))
    end
    all(v -> all(isfinite, v), inputs) || throw(ArgumentError("initial inputs must be finite"))
    return _validate_factor(T, factor, Fx, buffers)
end

function _validate_factor(::Type{T}, factor::Real, Fx, buffers) where {T}
    typed_factor = T(factor)
    isfinite(typed_factor) && typed_factor < zero(T) ||
        throw(ArgumentError("factor must be finite and negative after conversion"))
    return typed_factor
end

function _validate_factor(::Type{T}, factor::AbstractVector, Fx, buffers) where {T}
    eltype(factor) === T || throw(ArgumentError("factor vector must have the configured element type $T"))
    axes(factor) == axes(Fx) || throw(DimensionMismatch("factor and objective axes must match"))
    for buffer in buffers
        Base.mightalias(buffer, factor) && throw(ArgumentError("trial buffers must not alias factor"))
    end
    all(v -> isfinite(v) && v < zero(T), factor) ||
        throw(ArgumentError("factor entries must be finite and negative"))
    return factor
end

_validate_factor(::Type, factor, Fx, buffers) =
    throw(ArgumentError("factor must be a real number or a vector"))

@inline _factor_at(factor::Real, i) = factor
@inline _factor_at(factor::AbstractVector, i) = factor[i]

_validate_reference(::Type, ::Nothing, Fx, buffers) = Fx

function _validate_reference(::Type{T}, reference::AbstractVector, Fx, buffers) where {T}
    eltype(reference) === T ||
        throw(ArgumentError("reference must have the configured element type $T"))
    axes(reference) == axes(Fx) || throw(DimensionMismatch("reference and objective axes must match"))
    for buffer in buffers
        Base.mightalias(buffer, reference) && throw(ArgumentError("trial buffers must not alias reference"))
    end
    for i in eachindex(Fx)
        isfinite(reference[i]) && reference[i] >= Fx[i] ||
            throw(ArgumentError("reference must be finite and at least Fx in every component"))
    end
    return reference
end

_validate_reference(::Type, reference, Fx, buffers) =
    throw(ArgumentError("reference must be nothing or a vector"))

_initial_step(::Type, options, ::Nothing) = options.initial_step.alpha

function _initial_step(::Type{T}, options, alpha0::Real) where {T}
    alpha = T(alpha0)
    isfinite(alpha) && options.minstep <= alpha <= options.maxstep ||
        throw(ArgumentError("alpha0 must be finite and in [minstep, maxstep] after conversion"))
    return alpha
end

_initial_step(::Type, options, alpha0) = throw(ArgumentError("alpha0 must be nothing or a real number"))

function _step_limit(::Type{T}, limit) where {T}
    limit isa Real && !isnan(limit) && limit >= 0 ||
        throw(ArgumentError("maxstep must be a nonnegative real number or +Inf"))
    isinf(limit) && return T(Inf)
    converted = T(limit)
    # Respect the original bound even if conversion rounds upward or overflows.
    return converted > limit ? prevfloat(converted) : converted
end

@inline function _finish(search, alpha, reason, message, nfev, iterations)
    if search.options.verbose == 1
        @info "BackTracking finished" reason alpha nfev iterations
    end
    return LineSearchResult(reason === :success, alpha, reason, message,
                            nfev, 0, iterations, search.x_trial, search.F_trial)
end

"""
    perform_linesearch!(search::PreparedBackTracking, x, d, Fx;
                        slopes, factor, maxstep=Inf, reference=nothing, alpha0=nothing)

Accept a step when every trial value is finite and at most
`R[i] + c_1*alpha*η_i`, where `R = reference` (default `Fx`) and `η_i` is
`factor` if it is a scalar or `factor[i]` if it is a vector; every `η_i` must be
finite and negative. A `reference` must be finite and at least `Fx` in every
component. `slopes` supplies one finite slope per objective for interpolation
only; interpolation uses `Fx`, not `reference`.

`alpha0` replaces the configured initial step for this call and must lie in
`[minstep, maxstep]`. `maxstep` further bounds the initial trial. A bound below
`minstep` returns `:step_limit` without evaluation. See
[BackTracking](@ref backtracking) for reduction and stopping rules.
"""
function perform_linesearch!(search::PreparedBackTracking{T},
                             x::AbstractVector, d::AbstractVector,
                             Fx::AbstractVector; slopes, factor, maxstep=Inf,
                             reference=nothing, alpha0=nothing) where {T}
    options = search.options
    x_trial, F_trial = search.x_trial, search.F_trial
    factor_T = _validate_inputs(T, x, d, Fx, slopes, x_trial, F_trial, factor)
    reference_T = _validate_reference(T, reference, Fx, (x_trial, F_trial))
    initial = _initial_step(T, options, alpha0)
    limit = _step_limit(T, maxstep)
    # Preparation and _initial_step guarantee minstep <= initial <= the fixed maxstep.
    # Only the consumer bound can further restrict this initial trial/domain.
    if limit < options.minstep
        return _finish(search, limit, :step_limit,
                       "The consumer step limit is below the minimum", 0, 0)
    end
    alpha = min(initial, limit)
    nfev = 0
    iterations = 0
    while true
        @. x_trial = x + alpha * d
        if x_trial == x
            return _finish(search, alpha, :indistinguishable_point,
                           "The trial point is indistinguishable from x", nfev, iterations)
        end
        nfev += 1
        message = search.evaluate!(F_trial, x_trial)
        if message !== nothing
            message isa AbstractString ||
                throw(ArgumentError("evaluate! must return nothing or an error message"))
            return _finish(search, alpha, :evaluation_failure, String(message), nfev, iterations)
        end

        failed = nothing
        scale = options.c_1 * alpha
        for i in eachindex(Fx)
            if !isfinite(F_trial[i]) || !(F_trial[i] <= reference_T[i] + scale * _factor_at(factor_T, i))
                failed = i
                break
            end
        end
        if failed === nothing
            return _finish(search, alpha, :success, nothing, nfev, iterations)
        elseif iterations == options.max_iters
            return _finish(search, alpha, :max_iterations,
                           "Maximum number of backtracking reductions reached", nfev, iterations)
        end

        if any(v -> !isfinite(v), F_trial)
            alpha *= T(0.5)
        else
            j = failed
            alpha = _quadratic_step(alpha, F_trial[j], slopes[j], Fx[j],
                                    options.tau_1, options.tau_2)
        end
        iterations += 1
        if alpha <= options.minstep
            return _finish(search, alpha, :minimum_step,
                           "The reduced step is at or below the minimum", nfev, iterations)
        end
    end
end
