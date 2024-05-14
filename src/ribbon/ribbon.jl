export ribbon, ribbon!

using Backboner
using Makie

@recipe(Ribbon, chains) do scene
    Attributes(
        backgroundcolor = :black,
        colormap = :jet,
        colors = nothing,
        camcontrols = (;),
    )
end

include("utils.jl")
include("secondarystructure.jl")
include("render/render.jl")

function Makie.plot!(ribbon::Ribbon{Vector{Protein.Chain}})
    chains = deepcopy(ribbon[1])
    _assign_secondary_structure!(chains)

    isnothing(ribbon.colors) && (ribbon.colors = [LinRange(0, 1, length(chain)) for chain in protein])

    render!(ribbon, chains)

    return ribbon
end



#="""
    ribbon(protein::AbstractVector{Protein.Chain}; backgroundcolor=:black, camcontrols=(;), kwargs...)

Render a protein as a ribbon diagram. The display will be automatically centered on `protein`,
unless the user supplies `camcontrols` (see Makie's camera documentation for details).

See `render!` for additional keyword arguments.
"""
function ribbon(protein::AbstractVector{Protein.Chain}; backgroundcolor=:black, camcontrols=(;), kwargs...)
    scene = Scene(backgroundcolor=backgroundcolor)
    cam3d!(scene; camcontrols...)
    ribbon!(scene, protein; kwargs...)
    if isempty(camcontrols)
        center!(scene)
    end
    display(scene)
    return scene
end=#

#ribbon!(container, pdb_file::String; kwargs...) = ribbon!(container, Protein.readpdb(pdb_file); kwargs...)