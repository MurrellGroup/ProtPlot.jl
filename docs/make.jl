using GLMakie
using ProtPlot
using Documenter, Literate

GLMakie.activate!(px_per_unit = 2)

LITERATE_INPUT = joinpath(@__DIR__, "..", "examples")
LITERATE_OUTPUT = OUTPUT = joinpath(@__DIR__, "src/generated")

for (root, _, files) ∈ walkdir(LITERATE_INPUT), file ∈ files
    endswith(".jl")(file) || continue
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, LITERATE_INPUT=>LITERATE_OUTPUT))[1]
    Literate.markdown(ipath, opath)
end

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
        "Plot types" => [
            "Ribbon" => "generated/ribbon.md",
            "AtomPlot" => "generated/atomplot.md",
            "SpatialGraphPlot" => "generated/spatialgraphplot.md",
            "Ramachandran" => "generated/ramachandran.md",
        ],
    ],
    doctest=true,
)

deploydocs(;
    repo="github.com/MurrellGroup/ProtPlot.jl",
    devbranch="main",
)