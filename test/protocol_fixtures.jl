module ProtocolFixtures

using MultiObjectiveLineSearches: AbstractLineSearchMethod, LineSearchResult
import MultiObjectiveLineSearches: configure_linesearch, perform_linesearch!

# Test-only protocol witnesses. Neither fixture is a production line search or
# an implementation of Wolfe, Hager–Zhang, or a published nonmonotone method.
Base.@kwdef mutable struct PreparationCounts
    preparations::Int = 0
    conversions::Int = 0
    validations::Int = 0
end

struct CountedReal{T<:Real} <: Real
    value::T
    conversions::Base.RefValue{Int}
end
Base.Float32(v::CountedReal) = (v.conversions[] += 1; Float32(v.value))
Base.Float64(v::CountedReal) = (v.conversions[] += 1; Float64(v.value))

Base.@kwdef struct DerivativeProbe{A<:Real} <: AbstractLineSearchMethod
    alpha::A = 1.0
    maxstep::A = 1.0
    minstep::A = 1e-15
    max_iters::Int = 8
end

Base.@kwdef struct HistoryProbe{A<:Real} <: AbstractLineSearchMethod
    alpha::A = 1.0
    maxstep::A = 1.0
    minstep::A = 1e-15
    max_iters::Int = 8
    capacity::Int = 3
end

struct TrialOptions{T}
    alpha::T
    maxstep::T
    minstep::T
    max_iters::Int
end

struct DerivativeState{O,V,D,J,X,F,G}
    options::O
    evaluate!::V
    derivative!::D
    value_derivative!::J
    x_trial::X
    F_trial::F
    D_base::G
    D_trial::G
end

mutable struct HistoryState{O,V,X,F,H}
    options::O
    evaluate!::V
    x_trial::X
    F_trial::F
    last_x::X
    last_F::F
    history::H
    head::Int
    used::Int
    accepted::Int
end

function setup(spec, evaluate!, x, Fx, counts)
    counts.preparations += 1
    T = eltype(x)
    T <: AbstractFloat && isconcretetype(T) && eltype(Fx) === T ||
        throw(ArgumentError("matching concrete floating element types required"))
    !isempty(x) && !isempty(Fx) && all(isfinite,x) && all(isfinite,Fx) ||
        throw(ArgumentError("finite nonempty prototypes required"))
    counts.conversions += 1
    alpha, hi, lo = T(spec.alpha), T(spec.maxstep), T(spec.minstep)
    counts.validations += 1
    all(isfinite,(alpha,hi,lo)) && zero(T) < lo < hi && lo <= alpha <= hi &&
        spec.max_iters >= 0 || throw(ArgumentError("invalid converted options"))
    y,F = similar(x),similar(Fx)
    axes(y)==axes(x) && axes(F)==axes(Fx) || throw(DimensionMismatch("buffer axes"))
    eltype(y)===T && eltype(F)===T || throw(ArgumentError("buffer element types"))
    for output in (y,F), input in (x,Fx)
        Base.mightalias(output,input) && throw(ArgumentError("aliased prototype"))
    end
    Base.mightalias(y,F) && throw(ArgumentError("aliased outputs"))
    applicable(evaluate!,F,y) || throw(ArgumentError("missing value capability"))
    return TrialOptions(alpha,hi,lo,spec.max_iters),y,F
end

function configure_linesearch(spec::DerivativeProbe, evaluate!, x::AbstractVector,
                              Fx::AbstractVector; derivative! = nothing,
                              value_derivative! = nothing, counts=PreparationCounts())
    derivative! === nothing && value_derivative! === nothing &&
        throw(ArgumentError("directional derivatives required at preparation"))
    options,y,F = setup(spec,evaluate!,x,Fx,counts)
    D0,D = similar(F),similar(F)
    derivative! === nothing || applicable(derivative!,D,y,x) ||
        throw(ArgumentError("invalid separate derivative capability"))
    value_derivative! === nothing || applicable(value_derivative!,F,D,y,x) ||
        throw(ArgumentError("invalid joint capability"))
    return DerivativeState(options,evaluate!,derivative!,value_derivative!,y,F,D0,D)
end

function configure_linesearch(spec::HistoryProbe, evaluate!, x::AbstractVector,
                              Fx::AbstractVector; derivative! = nothing,
                              value_derivative! = nothing, counts=PreparationCounts())
    options,y,F = setup(spec,evaluate!,x,Fx,counts)
    spec.capacity > 0 || throw(ArgumentError("positive history capacity required"))
    last_x,last_F = copy(x),copy(Fx)
    history = Matrix{eltype(Fx)}(undef,length(Fx),spec.capacity)
    for (row,i) in enumerate(eachindex(Fx))
        history[row,1] = Fx[i]
    end
    return HistoryState(options,evaluate!,y,F,last_x,last_F,history,1,1,0)
end

outputs(s::DerivativeState) = (s.x_trial,s.F_trial,s.D_base,s.D_trial)
outputs(s::HistoryState) = (s.x_trial,s.F_trial,s.last_x,s.last_F,s.history)

function check_inputs(s,x,d,Fx)
    T = eltype(s.x_trial)
    all(v -> eltype(v)===T,(x,d,Fx)) || throw(ArgumentError("changed element type"))
    !isempty(x) && !isempty(Fx) || throw(ArgumentError("empty inputs"))
    axes(x)==axes(d)==axes(s.x_trial) && axes(Fx)==axes(s.F_trial) ||
        throw(DimensionMismatch("changed axes"))
    all(v -> all(isfinite,v),(x,d,Fx)) || throw(ArgumentError("nonfinite input"))
    # Tuple foreach specializes each buffer type, including the history matrix.
    foreach(outputs(s)) do output
        for input in (x,d,Fx)
            Base.mightalias(output,input) && throw(ArgumentError("input aliases owned storage"))
        end
    end
    return nothing
end

function step_limit(T,bound)
    bound isa Real && !isnan(bound) && bound >= 0 || throw(ArgumentError("invalid bound"))
    isinf(bound) && return T(Inf)
    converted = T(bound)
    return converted > bound ? prevfloat(converted) : converted
end

function status(message)
    message === nothing && return nothing
    message isa AbstractString || throw(ArgumentError("callback must return nothing or string"))
    return String(message)
end

finish(s,alpha,reason,nf,nd,it; message=nothing) =
    LineSearchResult(reason===:success,alpha,reason,message,nf,nd,it,s.x_trial,s.F_trial)

function perform_linesearch!(s::DerivativeState,x::AbstractVector,d::AbstractVector,
                             Fx::AbstractVector; maxstep=Inf)
    check_inputs(s,x,d,Fx)
    # Actual base/direction containers can change since preparation.
    if s.derivative! !== nothing
        applicable(s.derivative!,s.D_base,x,d) && applicable(s.derivative!,s.D_trial,s.x_trial,d) ||
            throw(ArgumentError("derivative callback cannot accept new containers"))
    end
    if s.value_derivative! !== nothing
        applicable(s.value_derivative!,s.F_trial,s.D_base,x,d) &&
            applicable(s.value_derivative!,s.F_trial,s.D_trial,s.x_trial,d) ||
            throw(ArgumentError("joint callback cannot accept new containers"))
    end
    T = eltype(x)
    upper = min(s.options.maxstep,step_limit(T,maxstep))
    upper < s.options.minstep && return finish(s,upper,:step_limit,0,0,0)
    alpha = min(s.options.alpha,upper)
    nf=nd=it=0
    if s.derivative! === nothing
        nf+=1; nd+=1
        message = status(s.value_derivative!(s.F_trial,s.D_base,x,d))
    else
        nd+=1
        message = status(s.derivative!(s.D_base,x,d))
    end
    message === nothing || return finish(s,alpha,:evaluation_failure,nf,nd,it;message)
    all(isfinite,s.D_base) && (s.derivative! !== nothing || all(isfinite,s.F_trial)) ||
        return finish(s,alpha,:nonfinite_trial,nf,nd,it)
    all(v -> v < zero(T),s.D_base) || return finish(s,alpha,:not_descent,nf,nd,it)
    while true
        @. s.x_trial = x + alpha*d
        s.x_trial==x && return finish(s,alpha,:indistinguishable_point,nf,nd,it)
        if s.value_derivative! !== nothing
            nf+=1; nd+=1
            message = status(s.value_derivative!(s.F_trial,s.D_trial,s.x_trial,d))
        else
            nf+=1
            message = status(s.evaluate!(s.F_trial,s.x_trial))
            message === nothing || return finish(s,alpha,:evaluation_failure,nf,nd,it;message)
            all(isfinite,s.F_trial) || return finish(s,alpha,:nonfinite_trial,nf,nd,it)
            nd+=1
            message = status(s.derivative!(s.D_trial,s.x_trial,d))
        end
        message === nothing || return finish(s,alpha,:evaluation_failure,nf,nd,it;message)
        all(isfinite,s.F_trial) && all(isfinite,s.D_trial) ||
            return finish(s,alpha,:nonfinite_trial,nf,nd,it)
        # Explicit test predicate, no scalarization or Wolfe condition.
        if all(i -> s.F_trial[i] <= Fx[i] && s.D_trial[i] <= zero(T),eachindex(Fx))
            return finish(s,alpha,:success,nf,nd,it)
        end
        it==s.options.max_iters && return finish(s,alpha,:probe_iterations,nf,nd,it)
        alpha *= T(0.5)
        it+=1
        alpha <= s.options.minstep && return finish(s,alpha,:minimum_step,nf,nd,it)
    end
end

function perform_linesearch!(s::HistoryState,x::AbstractVector,d::AbstractVector,
                             Fx::AbstractVector; maxstep=Inf)
    check_inputs(s,x,d,Fx)
    x==s.last_x && Fx==s.last_F || throw(ArgumentError("base is not last committed state"))
    T = eltype(x)
    upper = min(s.options.maxstep,step_limit(T,maxstep))
    upper < s.options.minstep && return finish(s,upper,:step_limit,0,0,0)
    alpha = min(s.options.alpha,upper)
    nf=it=0
    while true
        @. s.x_trial = x + alpha*d
        s.x_trial==x && return finish(s,alpha,:indistinguishable_point,nf,0,it)
        nf+=1
        message = status(s.evaluate!(s.F_trial,s.x_trial))
        message === nothing || return finish(s,alpha,:evaluation_failure,nf,0,it;message)
        all(isfinite,s.F_trial) || return finish(s,alpha,:nonfinite_trial,nf,0,it)
        accepted = true
        for (row,i) in enumerate(eachindex(Fx))
            ceiling = s.history[row,1]
            for column in 2:s.used
                ceiling = max(ceiling,s.history[row,column])
            end
            accepted &= s.F_trial[i] <= ceiling
        end
        if accepted
            # No history allocation or reinsertion of the base on later calls.
            copyto!(s.last_x,s.x_trial)
            copyto!(s.last_F,s.F_trial)
            slot = mod1(s.head+1,size(s.history,2))
            for (row,i) in enumerate(eachindex(Fx))
                s.history[row,slot] = s.F_trial[i]
            end
            s.head=slot
            s.used=min(s.used+1,size(s.history,2))
            s.accepted+=1
            return finish(s,alpha,:success,nf,0,it)
        end
        it==s.options.max_iters && return finish(s,alpha,:probe_iterations,nf,0,it)
        alpha *= T(0.5)
        it+=1
        alpha <= s.options.minstep && return finish(s,alpha,:minimum_step,nf,0,it)
    end
end

# The only consumer loop: algorithm-independent, same function for all methods.
function consume!(search,x,Fx,steps; maxstep=1)
    results = []
    for (d,data) in steps
        result = perform_linesearch!(search,x,d,Fx; maxstep,data...)
        push!(results,result)
        result.success || break
        copyto!(x,result.x_trial)
        copyto!(Fx,result.F_trial)
    end
    return results
end

end # module
