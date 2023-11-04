using Documenter, EsempioJulia

makedocs(
    modules = [EsempioJulia],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "gianluca.mastrantonio",
    sitename = "EsempioJulia.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/GianlucaMastrantonio/EsempioJulia.jl.git",
    push_preview = true
)
