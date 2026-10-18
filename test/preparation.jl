@testset "Prepared backtracking ($T)" for T in (Float32, Float64)
    calls = Ref(0)
    evaluator = (F,y) -> (calls[] += 1; fill!(F,-y[1]); nothing)
    x, Fx, d, slopes = T[0], T[0,0], T[1], T[-1,-1]
    spec = BackTracking()
    unavailable = (args...) -> error("unused derivative callback")
    search = @inferred configure_linesearch(spec,evaluator,x,Fx;
        derivative! = unavailable, value_derivative! = unavailable)
    @test calls[] == 0
    other = configure_linesearch(spec,evaluator,x,Fx)
    @test search.x_trial !== other.x_trial && search.F_trial !== other.F_trial
    @test !Base.mightalias(search.x_trial,x) && !Base.mightalias(search.F_trial,Fx)
    r = perform_linesearch!(search,x,d,Fx;slopes,factor=-one(T))
    @test r.success && r.nfev == calls[] == 1 && r.ndev == 0
    @test r.x_trial === search.x_trial && r.F_trial === search.F_trial
    copyto!(x,r.x_trial)
    copyto!(Fx,r.F_trial)
    @test !applicable(configure_linesearch,spec,T)
    @test !applicable(perform_linesearch!,search,evaluator,x,d,Fx,slopes,
                     similar(x),similar(Fx))
    for (px,pF) in ((T[],Fx),(x,T[]),(T[NaN],Fx),(x,T[Inf,0]),
                    (Int[0],Int[0,0]),(x,zeros(Int,2)))
        @test_throws ArgumentError configure_linesearch(spec,evaluator,px,pF)
    end
    @test_throws ArgumentError configure_linesearch(spec,nothing,x,Fx)
    @test_throws ArgumentError configure_linesearch(spec,identity,x,Fx)
    @test_throws MethodError perform_linesearch!(search,x,d,Fx;
        slopes,factor=-one(T),unsupported=1)
end

@testset "Consumer step domain ($T)" for T in (Float32, Float64)
    calls = Ref(0)
    trials = T[]
    evaluator = (F,y) -> (calls[] += 1; push!(trials,y[1]); fill!(F,-y[1]); nothing)
    x, Fx, d, slopes = T[0], T[0,0], T[1], T[-1,-1]
    search = configure_linesearch(BackTracking(minstep=0.25),evaluator,x,Fx)
    for bound in (0, 0.125, BigFloat(nextfloat(zero(T)))/4)
        before = calls[]
        r = perform_linesearch!(search,x,d,Fx;slopes,factor=-one(T),maxstep=bound)
        @test !r.success && r.reason === :step_limit
        @test r.nfev == r.ndev == r.iterations == 0
        @test r.alpha <= bound && calls[] == before
    end
    for bound in (0.25, 0.5, Inf, big(10)^1000)
        # Independent cases: discard each successful preparation before restarting.
        search = configure_linesearch(BackTracking(minstep=0.25),evaluator,x,Fx)
        r = perform_linesearch!(search,x,d,Fx;slopes,factor=-one(T),maxstep=bound)
        @test r.success && r.alpha == min(T(bound),one(T))
        @test trials[end] == r.alpha <= bound
    end
    # Nearest conversion rounds upward: the original real bound is authoritative.
    bound = BigFloat(T(0.5)) - BigFloat(eps(T(0.5)))/4
    search = configure_linesearch(BackTracking(minstep=0.25),evaluator,x,Fx)
    r = perform_linesearch!(search,x,d,Fx;slopes,factor=-one(T),maxstep=bound)
    @test r.success && r.alpha === prevfloat(T(0.5)) && r.alpha <= bound
    search = configure_linesearch(BackTracking(minstep=0.25),evaluator,x,Fx)
    for bound in (-1, -Inf, NaN, "invalid")
        before = calls[]
        @test_throws ArgumentError perform_linesearch!(search,x,d,Fx;
            slopes,factor=-one(T),maxstep=bound)
        @test calls[] == before
    end
    bounded = configure_linesearch(BackTracking(initial_step=InitialStatic(alpha=0.75),
        maxstep=0.75),evaluator,x,Fx)
    r = perform_linesearch!(bounded,x,d,Fx;slopes,factor=-one(T),maxstep=0.9)
    @test r.success && r.alpha === T(0.75)
    rejected = configure_linesearch(BackTracking(max_iters=1),
        (F,y)->(fill!(F,T(2));nothing),x,Fx)
    r = perform_linesearch!(rejected,x,d,Fx;slopes,factor=-one(T),maxstep=0.25)
    @test !r.success && r.reason === :max_iterations && r.alpha <= T(0.25)
end

@testset "BigFloat prepared execution" begin
    setprecision(256) do
        exact = big"0.12345678901234567890123456789"
        x, Fx, d, slopes = BigFloat[0], BigFloat[0,0], BigFloat[1], BigFloat[-1,-1]
        search = configure_linesearch(BackTracking(c_1=exact),DecreasingObjectives(),x,Fx)
        r = perform_linesearch!(search,x,d,Fx;slopes,factor=-big"1")
        @test search.options.c_1 == exact
        @test r.success && r.alpha == one(BigFloat)
        @test r.F_trial == BigFloat[-1,-1] && precision(r.alpha) == 256
    end
end
