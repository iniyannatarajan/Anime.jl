#push!(LOAD_PATH,"../src/")
#using Glob
using Anime
using Glob
using Documenter
using Literate

DocMeta.setdocmeta!(Anime, :DocTestSetup, :(using Anime); recursive=true)

# Setup Literate
OUTDIR = joinpath(@__DIR__, "src", "examples")
SOURCE_FILES = Glob.glob("*.jl", joinpath(dirname(pathof(Anime)), "..", "examples"))
println(SOURCE_FILES)
foreach(fn -> Literate.markdown(fn, OUTDIR, documenter=true), SOURCE_FILES)

EXAMPLES = [joinpath("examples", "createdataset.md"),
            joinpath("examples", "computecoherency.md"),
            joinpath("examples", "loadobservation.md"),
            joinpath("examples", "addmodels.md"),
            joinpath("examples", "pipeline.md")
           ]

makedocs(;
    modules=[Anime],
    authors="Iniyan Natarajan",
    repo="https://github.com/iniyannatarajan/Anime.jl/blob/{commit}{path}#{line}",
    sitename="Anime.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://iniyannatarajan.github.io/Anime.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "install.md",
        "Components" => "components.md",
        "Instrument Models" => "instrumentmodels.md",
        "Tutorial" => EXAMPLES,
        "Anime API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/iniyannatarajan/Anime.jl",
    devbranch="main",
    branch="gh-pages",
)
