export ribbon, ribbon!

using Backboner
using Makie

@recipe(Ribbon, chains) do scene
    Attributes(
        backgroundcolor = :black,
        colormap = :jet,
        colors = nothing,
        
        coil_radius = 0.2,
        coil_spline_quality = 20,
        coil_slice_quality = 20,

        helix_radius = 1.0,
        helix_width = 1.0,
        helix_thickness = 0.25,
        helix_spline_quality = 20,
        helix_slice_quality = 20,

        strand_width = 2.0,
        strand_thickness = 0.5,
        strand_spline_quality = 20,
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
    return scene
end=#

#ribbon!(container, pdb_file::String; kwargs...) = ribbon!(container, Protein.readpdb(pdb_file); kwargs...)