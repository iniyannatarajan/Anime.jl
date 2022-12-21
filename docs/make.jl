using Anime
using Documenter

DocMeta.setdocmeta!(Anime, :DocTestSetup, :(using Anime); recursive=true)

makedocs(;
    modules=[Anime],
    authors="Iniyan Natarajan <iniyannatarajan@gmail.com> and contributors",
    repo="https://github.com/iniyannatarajan/Anime.jl/blob/{commit}{path}#{line}",
    sitename="Anime.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://iniyannatarajan.github.io/Anime.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/iniyannatarajan/Anime.jl",
    devbranch="main",
)
