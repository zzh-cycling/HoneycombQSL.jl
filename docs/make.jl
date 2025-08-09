using HoneycombQSL
using Documenter

DocMeta.setdocmeta!(HoneycombQSL, :DocTestSetup, :(using HoneycombQSL); recursive=true)

makedocs(;
    modules=[HoneycombQSL],
    authors="Zhaohui Zhi",
    sitename="HoneycombQSL.jl",
    format=Documenter.HTML(;
        canonical="https://zzh-cycling.github.io/HoneycombQSL.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/zzh-cycling/HoneycombQSL.jl",
    devbranch="main",
)
