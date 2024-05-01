module Ribbon

export ribbon, ribbon!

using Backboner
using GLMakie
using Colors
using ColorSchemes
using LinearAlgebra

default_colorscheme = colorschemes[:jet]

include("utils.jl")
include("secondarystructure.jl")
include("shapes/shapes.jl")
include("segment.jl")
include("render.jl")

"""
    ribbon!(container, protein::AbstractVector{Protein.Chain}; kwargs...)

Renders a protein as a ribbon diagram onto a container.

Keyword arguments:
- `colorscheme::ColorScheme = ColorSchemes.jet`: The color scheme to use for the ribbon.
- `color_vectors::Vector{<:Vector{<:Union{AbstractFloat, RGB}}} = [LinRange(0, 1, length(chain)) for chain in protein]`:
    The color vectors to use for each chain. The length of each vector must match the length of the corresponding chain.
    The vectors must be either vectors of real number between 0 and 1, or vectors of RGB colors.
"""
function ribbon!(
    container,
    protein::AbstractVector{Protein.Chain};
    colorscheme::Union{ColorScheme, Symbol} = default_colorscheme,
    color_vectors::AbstractVector{<:AbstractVector{<:Union{Real, RGB}}} = [LinRange(0, 1, length(chain)) for chain in protein],
    kwargs...
)
    protein_assigned = ASS.assign_secondary_structure!(deepcopy(protein))

    colorscheme isa Symbol && (colorscheme = colorschemes[colorscheme])
    if eltype(eltype(color_vectors)) <: Real
        color_vectors = [colorscheme[color_vector] for color_vector in color_vectors]
    end

    render!(container, protein_assigned; colorscheme=colorscheme, color_vectors=color_vectors, kwargs...)

    return container
end

"""
    ribbon(protein::AbstractVector{Protein.Chain}; kwargs...)

Renders a protein as a ribbon diagram.
See `render!` for keyword arguments.
"""
function ribbon(protein::AbstractVector{Protein.Chain}; backgroundcolor=:black, kwargs...)
    scene = Scene(backgroundcolor=backgroundcolor)
    cam3d!(scene)
    ribbon!(scene, protein; kwargs...)
    center!(scene)
    display(scene)
    return scene
end

ribbon!(container, pdb_file::String; kwargs...) = ribbon!(container, Protein.readpdb(pdb_file); kwargs...)
ribbon(pdb_file::String; kwargs...) = ribbon(Protein.readpdb(pdb_file); kwargs...)

end