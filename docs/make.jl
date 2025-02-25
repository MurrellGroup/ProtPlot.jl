using ProtPlot
using Documenter

DocMeta.setdocmeta!(
    ProtPlot,
    :DocTestSetup,
    :(using ProtPlot);
    recursive=true,
)

makedocs(;
    modules=[ProtPlot],
    authors="Anton Oresten <antonoresten@gmail.com> and contributors",
    sitename="ProtPlot.jl",
    format=Documenter.HTML(;
        canonical="https://MurrellGroup.github.io/ProtPlot.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    doctest=true,
)

deploydocs(;
    repo="github.com/MurrellGroup/ProtPlot.jl",
    devbranch="main",
)