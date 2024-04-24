using CoherentNoise
using Documenter

DocMeta.setdocmeta!(CoherentNoise, :DocTestSetup, :(using CoherentNoise); recursive=true)

makedocs(;
    modules=[CoherentNoise],
    authors="Yusheng Zhao <yushengzhao2020@outlook.com> and contributors",
    sitename="CoherentNoise.jl",
    format=Documenter.HTML(;
        canonical="https://exAClior.github.io/CoherentNoise.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/exAClior/CoherentNoise.jl",
    devbranch="main",
)
