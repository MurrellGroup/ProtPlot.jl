module Ribbon

using Backboner
using GLMakie
using Colors
using ColorSchemes
using LinearAlgebra

include("utils.jl")
include("secondarystructure.jl")
include("shapes/shapes.jl")
include("segment.jl")
include("render.jl")

"""
    ribbon!(container, protein::Protein; kwargs...)

Renders a protein as a ribbon diagram onto a container.

Keyword arguments:
- `colorscheme::ColorScheme = ColorSchemes.jet`: The color scheme to use for the ribbon.
- `color_vectors::Vector{<:Vector{<:Union{AbstractFloat, RGB}}} = [LinRange(0, 1, length(chain)) for chain in protein]`:
    The color vectors to use for each chain. The length of each vector must match the length of the corresponding chain.
    The vectors must be either vectors of real number between 0 and 1, or vectors of RGB colors.
"""
function ribbon!(container, protein::Protein; kwargs...)
    protein_assigned = if has_assigned_ss(protein) # the function should probably be called `fully_assigned_ss` or something -- "has" is ambiguous
        protein
    else
        @debug """The given protein has one or more residues with unassigned secondary structure.
        Copying the protein and assigning secondary structure to the copy...
        This can be avoided by calling `assign_secondary_structure!` on the protein or by modifying the ssvector field of the protein chains manually."""
        assign_secondary_structure(protein)
    end
    render!(container, protein_assigned; kwargs...)

    return container
end

"""
    ribbon(protein::Protein; kwargs...)

Renders a protein as a ribbon diagram.
See `render!` for keyword arguments.
"""
function ribbon(protein::Protein; backgroundcolor=:black, kwargs...)
    scene = Scene(backgroundcolor=backgroundcolor)
    cam3d!(scene)
    ribbon!(scene, protein; kwargs...)
    center!(scene)
    display(scene)
    return scene
end

end