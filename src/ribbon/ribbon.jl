export Ribbon, ribbon, ribbon!, ribbon_scene

using Backboner
using Makie

@recipe(Ribbon, chains) do scene
    Attributes(
        backgroundcolor = :black,
        colormap = :jet,
        colors = nothing,

        coil_diameter = 0.4,
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
        strand_arrow_head_length = 5,
        strand_arrow_head_width = 3.5,
    )
end

include("utils.jl")
include("secondarystructure.jl")
include("render.jl")

# TODO: observe chains and re-render when they change
function Makie.plot!(ribbon::Ribbon{Tuple{Vector{Protein.Chain}}})
    chains = ribbon[1]

    chains = deepcopy(chains[])
    _assign_secondary_structure!(chains)
    isnothing(ribbon.colors[]) && (ribbon.colors = [LinRange(0, 1, length(chain)) for chain in chains])
    render!(ribbon, chains)

    return ribbon
end

Makie.convert_arguments(::Type{<:Ribbon}, chain::Protein.Chain) = ([chain],)

Makie.convert_arguments(::Type{<:Ribbon}, pdbfile::AbstractString) = (readpdb(pdbfile),)

#="""
    ribbon(protein::AbstractVector{Protein.Chain}; backgroundcolor=:black, camcontrols=(;), kwargs...)

Render a protein as a ribbon diagram. The display will be automatically centered on `protein`,
unless the user supplies `camcontrols` (see Makie's camera documentation for details).

See `render!` for additional keyword arguments.
"""=#
function ribbon_scene(args...; backgroundcolor=:black, camcontrols=(;), kwargs...)
    scene = Scene(backgroundcolor=backgroundcolor)
    cam3d!(scene; camcontrols...)
    ribbon!(scene, args...; kwargs...)
    if isempty(camcontrols)
        center!(scene)
    end
    return scene
end
