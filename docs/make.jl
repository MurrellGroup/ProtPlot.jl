using WGLMakie

using ProtPlot
using Documenter, Literate

LITERATE_INPUT = joinpath(@__DIR__, "..", "examples")
LITERATE_OUTPUT = OUTPUT = joinpath(@__DIR__, "src/generated")

for (root, _, files) ∈ walkdir(LITERATE_INPUT), file ∈ files
    @show root, files
    # ignore non julia files
    splitext(file)[2] == ".jl" || continue
    # full path to a literate script
    ipath = joinpath(root, file)
    # generated output path
    opath = splitdir(replace(ipath, LITERATE_INPUT=>LITERATE_OUTPUT))[1]
    # generate the markdown file calling Literate
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
            "AtomPlot" => "generated/atomplot.md",
            "SpatialGraphPlot" => "generated/spatialgraphplot.md",
        ],
    ],
    doctest=true,
)

deploydocs(;
    repo="github.com/MurrellGroup/ProtPlot.jl",
    devbranch="main",
)