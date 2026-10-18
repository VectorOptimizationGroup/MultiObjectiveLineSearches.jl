using Documenter
using DocumenterCitations
using MultiObjectiveLineSearches

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style = :numeric,
)

makedocs(
    sitename = "MultiObjectiveLineSearches.jl",
    modules = [MultiObjectiveLineSearches],
    plugins = [bib],
    checkdocs = :exports,
    format = Documenter.HTML(
        prettyurls = true,
        edit_link = nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => [
            "Installation" => "installation.md",
            "Quick Start" => "quickstart.md",
        ],
        "Guides" => [
            "Using a line search" => "guides/usage.md",
            "Extending the protocol" => "guides/extending.md",
        ],
        "Line Searches" => [
            "BackTracking" => "backtracking.md",
        ],
        "API Reference" => "api.md",
        "References" => "references.md",
    ],
)

deploydocs(
    repo = "github.com/VectorOptimizationGroup/MultiObjectiveLineSearches.jl.git",
    devbranch = "main",
)
