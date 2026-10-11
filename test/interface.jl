@testset "Result and independence" begin
    result = LineSearchResult(true, 0.5f0, :success, nothing, 2, 1)
    @test result.alpha === 0.5f0
    @test result.nfev == 2 && result.iterations == 1
    @test !ismutabletype(typeof(result))
    @test !hasfield(typeof(result), :F_final)
    # Loading configuration does not load solver/backend packages.
    forbidden = ("MOSolvers", "JuMP", "Ipopt", "HiGHS", "StochasticSolvers")
    @test all(id.name ∉ forbidden for id in keys(Base.loaded_modules))
end
