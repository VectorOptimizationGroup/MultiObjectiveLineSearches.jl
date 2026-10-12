"""
    BackTracking(; initial_step=InitialStatic(), c_1=1e-4, max_iters=1000,
                 minstep=1e-15, maxstep=1.0, verbose=0)

Backtracking line search for an Armijo-type sufficient decrease condition,
with safeguarded quadratic interpolation.

[`configure_linesearch`](@ref) checks the options after conversion:
`0 < c_1 < 1`, finite `0 < minstep < maxstep`, an initial step in
`[minstep, maxstep]`, `max_iters >= 0` and `verbose` equal to 0 or 1.
`max_iters` limits the number of step reductions; zero allows only the initial
trial. Set `verbose=1` to log the final diagnostic with `@info`.

See [BackTracking](@ref backtracking) for the condition, its sources and the
reduction rules.
"""
struct BackTracking{I<:InitialStatic,C<:Real,L<:Real,U<:Real} <: AbstractLineSearchMethod
    initial_step::I
    c_1::C
    max_iters::Int
    minstep::L
    maxstep::U
    verbose::Int
end

function BackTracking(; initial_step::InitialStatic=InitialStatic(),
                      c_1::Real=1e-4, max_iters::Integer=1000,
                      minstep::Real=1e-15, maxstep::Real=1.0, verbose::Integer=0)
    return BackTracking(initial_step, c_1, Int(max_iters), minstep, maxstep, Int(verbose))
end

struct BacktrackingParams{T<:AbstractFloat}
    initial_step::InitialStatic{T}
    c_1::T
    max_iters::Int
    minstep::T
    maxstep::T
    verbose::Int
end

function configure_linesearch(spec::BackTracking, ::Type{T}) where {T<:AbstractFloat}
    isconcretetype(T) || throw(ArgumentError("execution type must be concrete"))
    spec.max_iters >= 0 || throw(ArgumentError("max_iters must be nonnegative"))
    spec.verbose in (0, 1) || throw(ArgumentError("verbose must be 0 or 1"))
    c_1, lo, hi, alpha = T(spec.c_1), T(spec.minstep), T(spec.maxstep), T(spec.initial_step.alpha)
    isfinite(c_1) && zero(T) < c_1 < one(T) ||
        throw(ArgumentError("c_1 must be finite and in (0, 1) after conversion"))
    all(isfinite, (lo, hi, alpha)) && zero(T) < lo < hi && lo <= alpha <= hi ||
        throw(ArgumentError("require finite 0 < minstep < maxstep and minstep <= alpha <= maxstep after conversion"))
    return BacktrackingParams(InitialStatic(alpha), c_1, spec.max_iters, lo, hi, spec.verbose)
end

"""
    _quadratic_step(alpha, value, slope, value0)

Propose a step from the quadratic matching `value0`, `slope` and `value` at
`alpha`. Use the proposal only when `slope < 0`, the denominator is finite and
nonzero, and the proposal is finite and in `[0.05*alpha, 0.95*alpha]`;
otherwise return `alpha/2`.
See [Quadratic interpolation and safeguards](@ref quadratic-interpolation).
"""
function _quadratic_step(alpha::T, value::T, slope::T, value0::T) where {T}
    if slope < zero(T)
        den = value - value0 - alpha * slope
        if !iszero(den) && isfinite(den)
            candidate = (-slope * alpha^2) / (T(2) * den)
            if isfinite(candidate) && T(0.05) * alpha <= candidate <= T(0.95) * alpha
                return candidate
            end
        end
    end
    return alpha * T(0.5)
end

function _validate_inputs(::Type{T}, x, d, Fx, slopes, x_trial, F_trial, factor) where {T}
    inputs = (x, d, Fx, slopes)
    buffers = (x_trial, F_trial)
    all(v -> eltype(v) === T, (inputs..., buffers...)) ||
        throw(ArgumentError("all vectors must have the configured element type $T"))
    !isempty(x) && !isempty(Fx) || throw(ArgumentError("input vectors must be nonempty"))
    axes(x) == axes(d) == axes(x_trial) ||
        throw(DimensionMismatch("point, direction and trial-point axes must match"))
    axes(Fx) == axes(slopes) == axes(F_trial) ||
        throw(DimensionMismatch("objective, slope and trial-objective axes must match"))
    for buffer in buffers, input in inputs
        Base.mightalias(buffer, input) && throw(ArgumentError("trial buffers must not alias inputs"))
    end
    Base.mightalias(x_trial, F_trial) && throw(ArgumentError("trial buffers must not alias each other"))
    all(v -> all(isfinite, v), inputs) || throw(ArgumentError("initial inputs must be finite"))
    factor isa Real || throw(ArgumentError("factor must be real"))
    typed_factor = T(factor)
    isfinite(typed_factor) && typed_factor < zero(T) ||
        throw(ArgumentError("factor must be finite and negative after conversion"))
    return typed_factor
end

function _finish(search, alpha, reason, message, nfev, iterations)
    if search.verbose == 1
        @info "BackTracking finished" reason alpha nfev iterations
    end
    return LineSearchResult(reason === :success, alpha, reason, message, nfev, iterations)
end

"""
    perform_linesearch!(search::BacktrackingParams, evaluate!, x, d, Fx, slopes,
                        x_trial, F_trial; factor)

Search for a step `alpha` such that, for every objective `i`, the value at
`x + alpha*d` is finite and at most `Fx[i] + c_1*alpha*factor`. `factor` must
be finite and negative. `slopes` holds one slope per objective and is used only
to propose reduced steps.

If a trial has a nonfinite value, the step is halved; otherwise a rejected step
is reduced by safeguarded quadratic interpolation. The search fails when the
evaluator reports an error, after `max_iters` reductions, when a reduced step
is at or below `minstep`, or when `x + alpha*d` equals `x`. Evaluator
exceptions propagate.
See [BackTracking](@ref backtracking) for the acceptance and reduction rules.
"""
function perform_linesearch!(search::BacktrackingParams{T}, evaluate!,
                             x::AbstractVector, d::AbstractVector,
                             Fx::AbstractVector, slopes::AbstractVector,
                             x_trial::AbstractVector, F_trial::AbstractVector;
                             factor) where {T}
    factor_T = _validate_inputs(T, x, d, Fx, slopes, x_trial, F_trial, factor)
    alpha = search.initial_step.alpha
    nfev = 0
    iterations = 0
    while true
        @. x_trial = x + alpha * d
        if x_trial == x
            return _finish(search, alpha, :indistinguishable_point,
                           "The trial point is indistinguishable from x", nfev, iterations)
        end
        nfev += 1
        message = evaluate!(F_trial, x_trial)
        if message !== nothing
            message isa AbstractString ||
                throw(ArgumentError("evaluate! must return nothing or an error message"))
            return _finish(search, alpha, :evaluation_failure, String(message), nfev, iterations)
        end

        failed = nothing
        decrease = search.c_1 * alpha * factor_T
        for i in eachindex(Fx)
            if !isfinite(F_trial[i]) || !(F_trial[i] <= Fx[i] + decrease)
                failed = i
                break
            end
        end
        if failed === nothing
            return _finish(search, alpha, :success, nothing, nfev, iterations)
        elseif iterations == search.max_iters
            return _finish(search, alpha, :max_iterations,
                           "Maximum number of backtracking reductions reached", nfev, iterations)
        end

        if any(v -> !isfinite(v), F_trial)
            alpha *= T(0.5)
        else
            j = failed
            alpha = _quadratic_step(alpha, F_trial[j], slopes[j], Fx[j])
        end
        iterations += 1
        if alpha <= search.minstep
            return _finish(search, alpha, :minimum_step,
                           "The reduced step is at or below the minimum", nfev, iterations)
        end
    end
end
