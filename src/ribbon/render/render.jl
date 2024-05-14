include("shapes/shapes.jl")

function get_surface_vertices(ribbon::Ribbon, ::Val{:Coil}, segment_range::UnitRange{Int}, chain::Protein.Chain)
    adjusted_start = segment.range.start == 1 ? 1 : segment.range.start - 1
    adjusted_end = segment.range.stop == length(chain) ? length(chain) : segment.range.stop + 1
    points = Protein.alphacarbon_coords(chain)[:, adjusted_start:adjusted_end]
    attr = ribbon.attributes
    return coil_surface(points;
        attr.coil_radius, attr.coil_spline_quality, attr.coil_slice_quality)
end

function get_surface_vertices(ribbon::Ribbon, ::Val{:Helix}, segment_range::UnitRange{Int}, chain::Protein.Chain)
    points = Protein.alphacarbon_coords(chain)[:, segment_range] # startpoint is first point instead. including first N *and* CA could mess with normals
    return helix_surface(points; ribbon.attributes)
end

function get_surface_vertices(ribbon::Ribbon, ::Val{:Strand}, segment_range::UnitRange{Int}, chain::Protein.Chain)
    oxygen_coords = Protein.oxygen_coords(chain)[segment_range]
    oxygen_coords_side1 = oxygen_coords[:, begin:2:end-1]
    oxygen_coords_side2 = oxygen_coords[:, begin+1:2:end]
    points1 = hcat(startpoint, oxygen_coords_side1, endpoint)
    points2 = hcat(startpoint, oxygen_coords_side2, endpoint)
    attr = ribbon.attributes
    return arrow_surface(points1, points2;
        attr.strand_width, attr.strand_thickness, attr.strand_arrow_head_length, attr.strand_arrow_head_width, attr.strand_spline_quality)
end

function render!(ribbon::Ribbon, chain::Protein.Chain, colors::AbstractVector{<:RGB})
    for (ss_name, segment_range) in segments(chain)
        surface_vertices = get_surface_vertices(ribbon, Val(ss_name), segment_range, chain)
        color_matrix = expand_colors(colors, size(surface_vertices, 2))
        surface!(ribbon, eachslice(surface_vertices, dims=1)..., color=color_matrix)
    end
    return ribbon
end

function draw_lines_between_subchains!(ribbon::Ribbon, subchains::AbstractVector{Protein.Chain}, color::RGB; linewidth=2, kwargs...)
    for (i, j) in zip(1:length(subchains)-1, 2:length(subchains))
        startpoint, endpoint = subchains[i].backbone[end], subchains[j].backbone[begin]
        n_segments = trunc(Int, norm(endpoint - startpoint) / 0.8)
        xs, ys, zs = [LinRange(startpoint[i], endpoint[i], 2*n_segments) for i in 1:3]
        linesegments!(ribbon, xs, ys, zs; linewidth=linewidth, color=color, transparency=true)
    end
    return ribbon
end

function render!(ribbon::Ribbon, chains::Vector{Protein.Chain})
    for (chain, colors) in zip(chains, ribbon.attributes.color_vectors)
        subchain_ranges = split(chain)
        if length(subchain_ranges) == 1
            render!(ribbon, chain, colors; kwargs...)
        else
            subchains = [chain[r] for r in subchain_ranges]
            for (subchain, r) in zip(subchains, subchain_ranges)
                render!(ribbon, subchain, colors[r]; kwargs...)
            end
            draw_lines_between_subchains!(ribbon, subchains, missing_residue_color; kwargs...)
        end
    end
    return ribbon
end