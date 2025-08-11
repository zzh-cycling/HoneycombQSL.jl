using LatticeQSL
using Documenter

DocMeta.setdocmeta!(LatticeQSL, :DocTestSetup, :(using LatticeQSL); recursive=true)

makedocs(;
    modules=[LatticeQSL],
    authors="Zhaohui Zhi",
    sitename="LatticeQSL.jl",
    format=Documenter.HTML(;
        canonical="https://zzh-cycling.github.io/LatticeQSL.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/zzh-cycling/LatticeQSL.jl",
    devbranch="main",
)
