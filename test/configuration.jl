@testset "Typed configuration ($T)" for T in (Float32, Float64)
    spec = BackTracking()
    search = @inferred configure_linesearch(spec, T)
    @test search.c_1 === T(1e-4)
    @test search.tau_1 === T(0.05) && search.tau_2 === T(0.95)
    @test search.initial_step.alpha === one(T)
    @test search.minstep === T(1e-15)
    @test search.maxstep === one(T)
    @test search.max_iters == 1000
    @test search.verbose == 0
    @test isconcretetype(typeof(search))
    @test !ismutabletype(typeof(spec)) && !ismutabletype(typeof(search))
    @test !ismutabletype(typeof(search.initial_step))
    @test configure_linesearch(BackTracking(max_iters=0), T).max_iters == 0
    @test configure_linesearch(BackTracking(initial_step=InitialStatic(alpha=0.25), minstep=0.25), T).initial_step.alpha === T(0.25)
    # Configuration reuse across types cannot change the original specification.
    configure_linesearch(spec, Float64)
    @test configure_linesearch(spec, T) == search
    @test spec.initial_step.alpha === 1.0
    for value in (0.0, -1.0, 1.0, 2.0, NaN, Inf, -Inf)
        @test_throws ArgumentError configure_linesearch(BackTracking(c_1=value), T)
    end
    for value in (0.0, -1.0, NaN, Inf, -Inf, 2.0)
        @test_throws ArgumentError configure_linesearch(BackTracking(initial_step=InitialStatic(alpha=value)), T)
    end
    for value in (0.0, -1.0, NaN, Inf, -Inf, 1.0)
        @test_throws ArgumentError configure_linesearch(BackTracking(minstep=value), T)
    end
    for value in (0.0, -1.0, NaN, Inf, -Inf, 0.5)
        @test_throws ArgumentError configure_linesearch(BackTracking(maxstep=value), T)
    end
    halving = configure_linesearch(BackTracking(tau_1=0.5, tau_2=0.5), T)
    @test halving.tau_1 === halving.tau_2 === T(0.5)
    for value in (0.0, -0.1, 0.6, NaN, Inf, -Inf)
        @test_throws ArgumentError configure_linesearch(BackTracking(tau_1=value), T)
    end
    for value in (0.4, 1.0, 1.5, NaN, Inf, -Inf)
        @test_throws ArgumentError configure_linesearch(BackTracking(tau_2=value), T)
    end
    @test_throws ArgumentError configure_linesearch(BackTracking(max_iters=-1), T)
    @test_throws ArgumentError configure_linesearch(BackTracking(verbose=2), T)
end

@testset "Conversion boundaries" begin
    @test_throws ArgumentError configure_linesearch(BackTracking(c_1=1e-100), Float32)
    @test_throws ArgumentError configure_linesearch(BackTracking(c_1=prevfloat(1.0)), Float32)
    @test_throws ArgumentError configure_linesearch(BackTracking(tau_1=1e-100), Float32)
    @test_throws ArgumentError configure_linesearch(BackTracking(tau_2=prevfloat(1.0)), Float32)
    @test_throws ArgumentError configure_linesearch(BackTracking(minstep=1e-100), Float32)
    @test_throws ArgumentError configure_linesearch(BackTracking(maxstep=1e100), Float32)
    @test_throws ArgumentError configure_linesearch(BackTracking(minstep=1.0, maxstep=nextfloat(1.0)), Float32)
    @test_throws ArgumentError configure_linesearch(BackTracking(), AbstractFloat)
    exact = big"0.12345678901234567890123456789"
    @test BackTracking(c_1=exact).c_1 == exact
    @test configure_linesearch(BackTracking(c_1=exact), BigFloat).c_1 == exact
end
