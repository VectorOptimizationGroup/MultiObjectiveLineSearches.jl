@testset "Result and independence" begin
    result = LineSearchResult(true, 0.5f0, :success, nothing, 2, 0, 1, Float32[0.5], Float32[-0.5])
    @test result.alpha === 0.5f0
    @test result.nfev == 2 && result.iterations == 1
    @test !ismutabletype(typeof(result))
    @test result.ndev == 0
    @test !hasfield(typeof(result), :njoint)
    @test result.x_trial == Float32[0.5] && result.F_trial == Float32[-0.5]
    @test !hasfield(typeof(result), :F_final)
    # The package depends only on Julia's standard library.
    project = TOML.parsefile(joinpath(pkgdir(MultiObjectiveLineSearches), "Project.toml"))
    @test all(isdir(joinpath(Sys.STDLIB, name)) for name in keys(get(project, "deps", Dict())))
end
