@testset "Result and independence" begin
    result = LineSearchResult(true, 0.5f0, :success, nothing, 2, 1)
    @test result.alpha === 0.5f0
    @test result.nfev == 2 && result.iterations == 1
    @test !ismutabletype(typeof(result))
    @test !hasfield(typeof(result), :F_final)
    # The package depends only on Julia's standard library.
    project = TOML.parsefile(joinpath(pkgdir(MultiObjectiveLineSearches), "Project.toml"))
    @test all(isdir(joinpath(Sys.STDLIB, name)) for name in keys(get(project, "deps", Dict())))
end
