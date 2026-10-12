function search_case(T, evaluator; x=T[0], d=T[1], Fx=T[0, 0], slopes=T[-1, -1],
                     factor=-one(T), options=BackTracking())
    y, Fy = similar(x), similar(Fx)
    result = perform_linesearch!(configure_linesearch(options, T), evaluator,
                                 x, d, Fx, slopes, y, Fy; factor)
    return result, y, Fy
end

struct DecreasingObjectives end
function (::DecreasingObjectives)(F, x)
    fill!(F, -x[1])
    return nothing
end

@testset "Nonfinite trials and errors ($T)" for T in (Float32, Float64)
    for bad in T[NaN, Inf, -Inf], component in (1, 2), threshold in (T(0.75), T(0.3))
        trials = T[]
        function evaluate!(F, y)
            push!(trials, y[1])
            fill!(F, -y[1])
            y[1] > threshold && (F[component] = bad)
            return nothing
        end
        r, y, F = search_case(T, evaluate!)
        @test r.success && r.reason === :success
        @test r.alpha === y[1]
        @test F == fill(-r.alpha, 2)
        @test all(F .<= -T(1e-4)*r.alpha)
        @test trials == (threshold == T(0.75) ? T[1, 0.5] : T[1, 0.5, 0.25])
        @test r.nfev == length(trials) == r.iterations + 1
    end
    for bad in T[NaN, Inf, -Inf], component in (1, 2), limit in (0, 2)
        calls = Ref(0)
        function always_bad!(F, y)
            calls[] += 1
            fill!(F, -y[1]); F[component] = bad
            nothing
        end
        r, _, _ = search_case(T, always_bad!; options=BackTracking(max_iters=limit))
        @test !r.success && r.reason === :max_iterations
        @test r.nfev == calls[] == limit + 1
        @test r.iterations == limit
        @test r.alpha == T(0.5)^limit
    end
    trials = T[]
    function later_bad!(F, y)
        push!(trials, y[1])
        fill!(F, length(trials) == 1 ? one(T) : length(trials) == 2 ? T(NaN) : -y[1])
        nothing
    end
    r, _, F = search_case(T, later_bad!)
    @test r.success && r.nfev == 3 && r.iterations == 2
    @test trials[3] == trials[2]/2 && all(isfinite, F)
    for reduced in (false, true)
        evaluator = (F,y) -> reduced && y[1] == 1 ? (fill!(F,T(NaN)); nothing) : "reported error"
        r, _, _ = search_case(T, evaluator)
        @test !r.success && r.reason === :evaluation_failure && r.message == "reported error"
        @test r.nfev == (reduced ? 2 : 1)
        @test r.iterations == (reduced ? 1 : 0)
    end
    @test_throws ErrorException search_case(T, (F,y) -> error("unexpected"))
    @test_throws ArgumentError search_case(T, (F,y) -> 42)
end

@testset "Finite interpolation and limits ($T)" for T in (Float32, Float64)
    quadratic!(F,y) = (F[1] = (y[1]-one(T))^2; nothing)
    r,y,F = search_case(T, quadratic!; d=T[3], Fx=T[1], slopes=T[-6], factor=T(-6))
    @test r.success && r.nfev == 2 && r.iterations == 1
    @test r.alpha ≈ one(T)/T(3)
    @test y ≈ T[1]
    @test F ≈ T[0] atol=eps(T)
    @test r.alpha isa T && eltype(F) === T
    # The proposal 1/3 is used only inside [tau_1*alpha, tau_2*alpha]; otherwise halve.
    for (tau_1, tau_2, alpha) in ((0.3, 0.7, 1/3), (0.4, 0.95, 0.5), (0.5, 0.5, 0.5))
        r,_,_ = search_case(T, quadratic!; d=T[3], Fx=T[1], slopes=T[-6], factor=T(-6),
                            options=BackTracking(; tau_1, tau_2))
        @test r.success && r.nfev == 2 && r.alpha ≈ T(alpha)
    end
    # Both components violate; the first determines interpolation, independently of factor.
    trials = T[]
    function first_violation!(F,y)
        push!(trials,y[1]); fill!(F, length(trials)==1 ? T(4) : T(-10)); nothing
    end
    r,_,_ = search_case(T, first_violation!; Fx=T[1,1], slopes=T[-6,-12], factor=T(-2))
    @test r.success && trials ≈ T[1, 1/3]
    # Upper bound: the proposal 0.8 is kept with tau_2 = 0.95 and halved with tau_2 = 0.6.
    upper!(F,y) = (push!(trials,y[1]); F[1] = length(trials)==1 ? T(-0.375) : T(-10); nothing)
    for (tau_2, step) in ((0.95, 0.8), (0.6, 0.5))
        empty!(trials)
        r,_,_ = search_case(T, upper!; Fx=T[0], slopes=T[-1], factor=-one(T),
                            options=BackTracking(; c_1=0.5, tau_2))
        @test r.success && trials ≈ T[1, step]
    end
    # Nonnegative component slope is valid and uses halving, not scalar descent rejection.
    empty!(trials)
    r,_,_ = search_case(T, first_violation!; Fx=T[1,1], slopes=T[0, -12], factor=T(-2))
    @test r.success && trials == T[1,0.5]
    for limit in (0,2)
        r,_,_ = search_case(T,(F,y)->(fill!(F,T(2));nothing); options=BackTracking(max_iters=limit))
        @test !r.success && r.reason === :max_iterations && r.nfev == limit+1
    end
    r,_,_ = search_case(T,(F,y)->(fill!(F,T(NaN));nothing); options=BackTracking(minstep=0.25))
    @test r.reason === :minimum_step && r.alpha == T(0.25)
    @test r.nfev == 2 && r.iterations == 2
    # Initial alpha equal to the minimum is still evaluated.
    r,_,_ = search_case(T,DecreasingObjectives(); options=BackTracking(initial_step=InitialStatic(alpha=0.25),minstep=0.25))
    @test r.success && r.alpha == T(0.25) && r.nfev == 1
    calls = Ref(0)
    evaluator=(F,y)->(calls[]+=1; fill!(F,T(NaN)); nothing)
    r,_,_ = search_case(T,evaluator; x=T[1],d=T[-1],options=BackTracking(minstep=1e-30))
    @test r.reason === :indistinguishable_point
    @test 0 < r.nfev == calls[] < 1000
    @test r.iterations == r.nfev
    @test one(T)-r.alpha == one(T)
    r,_,_ = search_case(T,DecreasingObjectives();d=T[0])
    @test r.reason === :indistinguishable_point && r.nfev == 0
end

@testset "State, buffers and verbosity ($T)" for T in (Float32,Float64)
    x,d,Fx,slopes = T[0],T[1],T[0,0],T[-1,-1]
    original = map(copy,(x,d,Fx,slopes))
    y,F = similar(x),similar(Fx)
    search=configure_linesearch(BackTracking(),T)
    r=@inferred perform_linesearch!(search,DecreasingObjectives(),x,d,Fx,slopes,y,F;factor=-one(T))
    @test r.success && r.nfev == 1 && r.iterations == 0
    @test y == T[1] && F == T[-1,-1]
    @test (x,d,Fx,slopes)==original
    snapshot=copy(F)
    r2=perform_linesearch!(search,(F,y)->(fill!(F,T(NaN));"stop"),x,d,Fx,slopes,y,F;factor=-one(T))
    @test !r2.success && (x,d,Fx,slopes)==original
    @test snapshot == T[-1,-1] && r.alpha == one(T)
    r3=perform_linesearch!(search,DecreasingObjectives(),x,d,Fx,slopes,y,F;factor=-one(T))
    @test r3 == r && F == snapshot
    noisy=configure_linesearch(BackTracking(verbose=1),T)
    @test_logs (:info,"BackTracking finished") begin
        r4=perform_linesearch!(noisy,DecreasingObjectives(),x,d,Fx,slopes,y,F;factor=-one(T))
        @test r4 == r
    end
end

@testset "Invalid inputs before evaluation ($T)" for T in (Float32,Float64)
    calls=Ref(0)
    evaluator=(F,y)->(calls[]+=1;fill!(F,zero(T));nothing)
    search=configure_linesearch(BackTracking(),T)
    function check(args...;factor=-one(T), error=ArgumentError)
        @test_throws error perform_linesearch!(search,evaluator,args...;factor)
        @test calls[] == 0
    end
    x,d,F,s = T[0],T[1],T[0,0],T[-1,-1]
    y,Fy=similar(x),similar(F)
    for factor in (0,1,Inf,-Inf,NaN,"invalid")
        check(x,d,F,s,y,Fy;factor)
    end
    check(T[],T[],F,s,T[],Fy)
    check(x,d,T[],T[],y,T[])
    check(x,T[1,2],F,s,y,Fy;error=DimensionMismatch)
    check(x,d,F,T[-1],y,Fy;error=DimensionMismatch)
    check(x,d,F,s,T[0,0],Fy;error=DimensionMismatch)
    check(x,d,F,s,y,T[0];error=DimensionMismatch)
    check(x,d,F,s,zeros(Int,1),Fy)
    for k in 1:4
        arrays=map(copy,(x,d,F,s)); arrays[k][1]=T(NaN)
        check(arrays...,y,Fy)
    end
    # Use equal lengths to isolate each alias rejection from dimensional checks.
    inputs=[T[0,0],T[1,1],T[0,0],T[-1,-1]]
    for input in inputs
        check(inputs...,view(input,:),zeros(T,2))
        check(inputs...,zeros(T,2),view(input,:))
    end
    buffer=zeros(T,2)
    check(inputs...,buffer,view(buffer,:))
    T===Float32 && check(x,d,F,s,y,Fy;factor=-1e-100)
end
