using Documenter
using DocumenterCitations
using MultiObjectiveLineSearches

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style = :numeric,
)

# No remote yet: `remotes = nothing` and no `deploydocs` until publication.
makedocs(
    sitename = "MultiObjectiveLineSearches.jl",
    modules = [MultiObjectiveLineSearches],
    plugins = [bib],
    checkdocs = :exports,
    remotes = nothing,
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
        ],
        "Line Searches" => [
            "BackTracking" => "backtracking.md",
        ],
        "API Reference" => "api.md",
        "References" => "references.md",
    ],
)
