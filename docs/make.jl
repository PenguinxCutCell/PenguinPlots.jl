using Documenter
using PenguinPlots

makedocs(
    modules=[PenguinPlots],
    authors="PenguinxCutCell contributors",
    sitename="PenguinPlots.jl",
    format=Documenter.HTML(
        canonical="https://PenguinxCutCell.github.io/PenguinPlots.jl",
        repolink="https://github.com/PenguinxCutCell/PenguinPlots.jl",
        collapselevel=2,
    ),
    doctest=false,
    pages=[
        "Home" => "index.md",
        "Scalar" => "scalar.md",
        "Flow" => "flow.md",
        "Stefan" => "stefan.md",
        "Coupling" => "coupling.md",
        "API" => "api.md",
    ],
    pagesonly=true,
    warnonly=false,
    remotes=nothing,
)

if get(ENV, "CI", "") == "true"
    deploydocs(
        repo="github.com/PenguinxCutCell/PenguinPlots.jl",
        push_preview=true,
    )
end
