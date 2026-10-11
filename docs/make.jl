using Documenter
using DocumenterCitations
using MultiObjectiveLineSearches

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

makedocs(
    sitename="MultiObjectiveLineSearches.jl",
    modules=[MultiObjectiveLineSearches],
    plugins=[bib],
    checkdocs=:exports,
    remotes=nothing,
    format=Documenter.HTML(prettyurls=false, edit_link=nothing),
    pages=[
        "Home" => "index.md",
        "Guides" => ["Using a line search" => "guides/usage.md"],
        "Line Searches" => ["BackTracking" => "backtracking.md"],
        "API Reference" => "api.md",
        "References" => "references.md",
    ],
)
