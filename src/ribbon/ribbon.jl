export Ribbon, ribbon, ribbon!, ribbon_scene

using Backboner
using Makie

@recipe(Ribbon, chains) do scene
    Attributes(
        secondary_structures = nothing,
        colors = nothing,
        colormap = :jet,

        coil_diameter = 0.4,
        coil_spline_quality = 20,
        coil_slice_quality = 20,

        helix_width = 2.0,
        helix_thickness = 0.5,
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
function Makie.plot!(ribbon::Ribbon{<:Tuple{<:AbstractVector{Protein.Chain}}})
    chains = ribbon[1][]
    isnothing(ribbon.secondary_structures[]) && (ribbon.secondary_structures[] = _assign_secondary_structure(chains))
    isnothing(ribbon.colors[]) && (ribbon.colors[] = [range(0, 1, length(chain)) for chain in chains])
    render!(ribbon, chains)
    return ribbon
end

Makie.convert_arguments(::Type{<:Ribbon}, chain::Protein.Chain) = ([chain],)

Makie.convert_arguments(::Type{<:Ribbon}, pdbfile::AbstractString) = (readpdb(pdbfile),)

Makie.convert_arguments(::Type{<:Ribbon}, backbones::AbstractVector{<:Backbone}) = (Protein.Chain.(backbones),)
Makie.convert_arguments(T::Type{<:Ribbon}, backbone::Backbone) = Makie.convert_arguments(T, [backbone])

Makie.convert_arguments(T::Type{<:Ribbon}, coords_vec::AbstractVector{<:AbstractMatrix{<:Real}}) = Makie.convert_arguments(T, Backbone.(coords_vec))
Makie.convert_arguments(T::Type{<:Ribbon}, coords::AbstractMatrix{<:Real}) = Makie.convert_arguments(T, [coords])

"""
    ribbon_scene(chains::AbstractVector{Protein.Chain}; backgroundcolor=:black, camcontrols=(;), kwargs...)

Render a protein as a ribbon diagram. The display will be automatically centered on the rendered ribbon,
unless the user supplies `camcontrols` (see Makie's camera documentation for details).
"""
function ribbon_scene(args...; backgroundcolor=:black, camcontrols=(;), kwargs...)
    scene = Scene(backgroundcolor=backgroundcolor)
    cam3d!(scene; camcontrols...)
    ribbon!(scene, args...; kwargs...)
    if isempty(camcontrols)
        center!(scene)
    end
    return scene
end
