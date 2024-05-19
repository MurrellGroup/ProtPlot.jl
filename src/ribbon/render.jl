include("shapes/shapes.jl")

function get_surface_vertices(ribbon::Ribbon, ::Val{:Coil}, segment_range::UnitRange{Int}, chain::Protein.Chain)
    adjusted_range = max(1, segment_range.start - 1):min(length(chain), segment_range.stop + 1)
    points = Protein.alphacarbon_coords(chain)
    return coil_surface(ribbon.attributes, points; segment_range=adjusted_range), adjusted_range
end

function get_surface_vertices(ribbon::Ribbon, ::Val{:Helix}, segment_range::UnitRange{Int}, chain::Protein.Chain)
    points = Protein.alphacarbon_coords(chain)
    return helix_surface(ribbon.attributes, points; segment_range), segment_range
end

function get_surface_vertices(ribbon::Ribbon, ::Val{:Strand}, segment_range::UnitRange{Int}, chain::Protein.Chain)
    ca_points = Protein.alphacarbon_coords(chain)
    o_points = Protein.oxygen_coords(chain) .|> eltype(ca_points)
    return strand_surface(ribbon.attributes, ca_points, o_points; segment_range), segment_range
end

function render!(ribbon::Ribbon, chain::Protein.Chain, color::AbstractVector{<:Real}, colormap)
    for (ss_name, segment_range) in segments(chain)
        surface_vertices, adjusted_range = get_surface_vertices(ribbon, Val(ss_name), segment_range, chain)
        surface!(ribbon, eachslice(surface_vertices, dims=1)...,
            color = expand_colors(color[adjusted_range], size(surface_vertices, 2)),
            colormap = colormap,
            colorrange = (0, 1))
    end
    return ribbon
end

function render!(ribbon::Ribbon, chain::Protein.Chain, color::AbstractVector{<:Colorant}, colormap)
    render!(ribbon, chain, LinRange(0, 1, length(color)), color)
end

function draw_lines_between_subchains!(
    ribbon::Ribbon, chain::Protein.Chain, subchain_ranges::Vector{UnitRange{Int}},
    color, linewidth = 2
)
    subchains = [chain[r] for r in subchain_ranges]
    for (i, j) in zip(1:length(subchains)-1, 2:length(subchains))
        startpoint, endpoint = subchains[i].backbone[end-1], subchains[j].backbone[begin+1]
        n_segments = trunc(Int, norm(endpoint - startpoint) / 0.8)
        xs, ys, zs = [LinRange(startpoint[i], endpoint[i], 2*n_segments) for i in 1:3]
        linesegments!(ribbon, xs, ys, zs; linewidth=linewidth, color=color, transparency=true)
    end
    return ribbon
end

function render!(ribbon::Ribbon, chains::Vector{Protein.Chain})
    for (chain, color) in zip(chains, ribbon.colors[])
        length(chain) > 1 || throw(ArgumentError("Chain must have at least 2 residues"))
        length(chain) == length(color) || throw(ArgumentError("Chain and color must have the same length"))
        subchain_ranges = get_subchain_ranges(chain)
        for r in subchain_ranges
            render!(ribbon, chain[r], color[r], ribbon.colormap)
        end
        draw_lines_between_subchains!(ribbon, chain, subchain_ranges, :lightgray)
    end
    return ribbon
end