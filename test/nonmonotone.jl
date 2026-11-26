@testset "Componentwise factor ($T)" for T in (Float32, Float64)
    # F_1 decreases fast and F_2 slowly along d: a common factor max(s) = -1 accepts
    # alpha = 1 only if both decrease by c_1; the componentwise term asks more of F_1.
    linear!(F, y) = (F[1] = -T(4) * y[1]; F[2] = -y[1]; nothing)
    x, d, Fx, slopes = T[0], T[1], T[0, 0], T[-4, -1]
    options = BackTracking(c_1=0.5)
    search = configure_linesearch(options, linear!, x, Fx)
    r = perform_linesearch!(search, x, d, Fx; slopes, factor=slopes)
    @test r.success && r.alpha == one(T) && r.nfev == 1
    # F_1 is now flat at the first trial: rejected by its own term, accepted by max(s).
    trials = T[]
    flat!(F, y) = (push!(trials, y[1]); F[1] = length(trials) == 1 ? -T(1) : -T(4) * y[1];
                   F[2] = -y[1]; nothing)
    search = configure_linesearch(options, flat!, x, Fx)
    r = perform_linesearch!(search, x, d, Fx; slopes, factor=maximum(slopes))
    @test r.success && trials == T[1]
    empty!(trials)
    r = perform_linesearch!(search, x, d, Fx; slopes, factor=slopes)
    @test r.success && r.nfev == 2 && trials[1] == 1 && trials[2] < 1
    # A vector with equal entries reproduces the scalar factor exactly.
    quadratic!(F, y) = (F[1] = (y[1] - 1)^2; F[2] = (y[1] - 1)^2 + 1; nothing)
    for tau in ((0.05, 0.95), (0.5, 0.5))
        bt = BackTracking(tau_1=tau[1], tau_2=tau[2])
        a = configure_linesearch(bt, quadratic!, T[0], T[1, 2])
        b = configure_linesearch(bt, quadratic!, T[0], T[1, 2])
        ra = perform_linesearch!(a, T[0], T[30], T[1, 2]; slopes=T[-60, -60], factor=T(-60))
        rb = perform_linesearch!(b, T[0], T[30], T[1, 2]; slopes=T[-60, -60], factor=T[-60, -60])
        @test (ra.success, ra.alpha, ra.nfev, ra.iterations) == (rb.success, rb.alpha, rb.nfev, rb.iterations)
        @test ra.x_trial == rb.x_trial && ra.F_trial == rb.F_trial
    end
end

@testset "Reference and per-call initial step ($T)" for T in (Float32, Float64)
    calls = Ref(0)
    up!(F, y) = (calls[] += 1; F[1] = y[1]; F[2] = -y[1]; nothing)   # F_1 increases along d
    x, d, Fx, slopes = T[0], T[1], T[0, 0], T[-1, -1]
    search = configure_linesearch(BackTracking(minstep=1e-3), up!, x, Fx)
    r = perform_linesearch!(search, x, d, Fx; slopes, factor=slopes)
    @test r.reason === :minimum_step
    r = perform_linesearch!(search, x, d, Fx; slopes, factor=slopes, reference=T[2, 0])
    @test r.success && r.alpha == one(T) && r.nfev == 1 && r.F_trial == T[1, -1]
    r = perform_linesearch!(search, x, d, Fx; slopes, factor=slopes, reference=Fx)
    @test r.reason === :minimum_step
    calls[] = 0
    for (reference, error) in ((T[-1, 0], ArgumentError), (T[NaN, 0], ArgumentError),
                               (T[Inf, 0], ArgumentError), (T[2], DimensionMismatch),
                               ([2.0f0, 0.0f0] .|> (T === Float32 ? Float64 : Float32), ArgumentError),
                               (search.F_trial, ArgumentError), (view(search.F_trial, :), ArgumentError),
                               (2, ArgumentError))
        @test_throws error perform_linesearch!(search, x, d, Fx; slopes, factor=slopes, reference)
    end
    @test calls[] == 0
    # Interpolation uses Fx, not the reference: the proposal stays 1/3.
    quadratic!(F, y) = (F[1] = (y[1] - 1)^2; nothing)
    q = configure_linesearch(BackTracking(), quadratic!, T[0], T[1])
    r = perform_linesearch!(q, T[0], T[3], T[1]; slopes=T[-6], factor=T[-6], reference=T[2])
    @test r.success && r.nfev == 2 && r.alpha ≈ one(T) / 3

    # Per-call initial step.
    wide = BackTracking(minstep=prevfloat(1e-20), maxstep=1e10)
    s = configure_linesearch(wide, DecreasingObjectives(), x, Fx)
    r = perform_linesearch!(s, x, d, Fx; slopes, factor=-one(T))
    @test r.success && r.alpha == one(T)
    r = perform_linesearch!(s, x, d, Fx; slopes, factor=-one(T), alpha0=1e5)
    @test r.success && r.alpha == T(1e5) && r.nfev == 1
    r = perform_linesearch!(s, x, d, Fx; slopes, factor=-one(T), alpha0=1e5, maxstep=4)
    @test r.success && r.alpha == T(4)
    r = perform_linesearch!(s, x, d, Fx; slopes, factor=-one(T), alpha0=1e-10)
    @test r.success && r.alpha == T(1e-10)
    calls[] = 0
    t = configure_linesearch(BackTracking(), up!, x, Fx)
    for alpha0 in (2, 1e-16, 0, -1, NaN, Inf, "1")
        @test_throws ArgumentError perform_linesearch!(t, x, d, Fx; slopes, factor=slopes, alpha0)
    end
    @test calls[] == 0
    for bad in (T[-1], T[-1, 0], T[-1, NaN], (T === Float32 ? [-1.0, -1.0] : Float32[-1, -1]),
                view(t.F_trial, :))
        @test_throws Union{ArgumentError,DimensionMismatch} perform_linesearch!(t, x, d, Fx; slopes, factor=bad)
    end
    @test calls[] == 0
end

# Literal transcription of the steepest descent loop of Fliege and Svaiter with the
# max-type nonmonotone test of Mita, Fukuda and Yamashita (Algorithm 3) and an initial
# step supplied by the caller; used to check that BackTracking reproduces it exactly.
module LiteralDescent

objectives(x) = [sum((x .- 1) .^ 2) + sum(sin.(3 .* x)) / 2,
                 sum(10 .* (x[2:end] .- x[1:end-1] .^ 2) .^ 2) + sum((1 .- x) .^ 2)]

function gradients(x)
    n = length(x)
    g1 = 2 .* (x .- 1) .+ 1.5 .* cos.(3 .* x)
    g2 = -2 .* (1 .- x)
    for i in 1:n-1
        r = x[i+1] - x[i]^2
        g2[i+1] += 20r
        g2[i] -= 40r * x[i]
    end
    return hcat(g1, g2)
end

# Minimum-norm element of the convex hull of the two columns of V.
function weights(V)
    g1, g2 = V[:, 1], V[:, 2]
    den = sum(abs2, g1 .- g2)
    w = den > 0 ? clamp(sum((g2 .- g1) .* g2) / den, 0.0, 1.0) : 0.5
    return [w, 1 - w]
end

function bb_step(x, xprev, V, Vprev, γ)
    isempty(xprev) && return 1.0
    s = x .- xprev
    y = (V .- Vprev) * γ
    sy = sum(s .* y)
    return sy > 0 ? clamp(sum(s .* s) / sy, 1e-10, 1e10) : 1.0
end

# One search: halve from t until the test holds or t < 1e-20.
function search(x, d, ref, slopes, t; c=1e-4)
    nf = 0
    red = 0
    while t >= 1e-20
        xt = x .+ t .* d
        Ft = objectives(xt)
        nf += 1
        all(Ft .<= ref .+ c * t .* slopes) && return (true, t, nf, red, xt, Ft)
        t /= 2
        red += 1
    end
    return (false, t, nf, red, x, objectives(x))
end

end # module

@testset "Exact agreement with a literal nonmonotone loop (memory $M, BB $bb)" for M in (0, 3), bb in (false, true)
    L = LiteralDescent
    x0 = [-1.2, 1.0, -0.5, 0.8]
    function run(use_package)
        x = copy(x0); Fx = L.objectives(x); history = [copy(Fx)]
        xprev = Float64[]; Vprev = zeros(0, 0)
        evaluate! = (F, y) -> (copyto!(F, L.objectives(y)); nothing)
        prepared = configure_linesearch(
            BackTracking(c_1=1e-4, tau_1=0.5, tau_2=0.5, minstep=prevfloat(1e-20),
                         maxstep=1e10, max_iters=10_000), evaluate!, x, Fx)
        trace = Tuple{Float64,Int,Int}[]
        for _ in 1:200            # stops earlier at a Pareto critical point
            V = L.gradients(x)
            γ = L.weights(V)
            d = -(V * γ)
            slopes = V' * d
            all(<(0), slopes) || break
            t0 = bb ? L.bb_step(x, xprev, V, Vprev, γ) : 1.0
            ref = reduce((a, b) -> max.(a, b), history)
            if use_package
                r = perform_linesearch!(prepared, x, d, Fx; slopes, factor=slopes,
                                        reference=ref, alpha0=t0)
                r.success || break
                push!(trace, (r.alpha, r.nfev, r.iterations))
                xnew, Fnew = copy(r.x_trial), copy(r.F_trial)
            else
                ok, t, nf, red, xnew, Fnew = L.search(x, d, ref, slopes, t0)
                ok || break
                push!(trace, (t, nf, red))
            end
            xprev, Vprev = x, V
            x, Fx = xnew, Fnew
            push!(history, copy(Fx))
            length(history) > M + 1 && popfirst!(history)
        end
        return x, Fx, trace
    end
    xa, Fa, ta = run(false)
    xb, Fb, tb = run(true)
    @test length(ta) >= 40
    @test ta == tb
    @test xa == xb && Fa == Fb
    @test sum(last, ta) > 0          # the comparison exercises reductions
end
