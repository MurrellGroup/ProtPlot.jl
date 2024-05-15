include("shapes/shapes.jl")

function get_surface_vertices(ribbon::Ribbon, ::Val{:Coil}, segment_range::UnitRange{Int}, chain::Protein.Chain)
    adjusted_start = segment_range.start == 1 ? 1 : segment_range.start - 1
    adjusted_end = segment_range.stop == length(chain) ? length(chain) : segment_range.stop + 1
    points = Protein.alphacarbon_coords(chain)
    return coil_surface(ribbon.attributes, points; segment_range=adjusted_start:adjusted_end)
end

function get_surface_vertices(ribbon::Ribbon, ::Val{:Helix}, segment_range::UnitRange{Int}, chain::Protein.Chain)
    points = Protein.alphacarbon_coords(chain)
    return helix_surface(ribbon.attributes, points; segment_range)
end

function get_surface_vertices(ribbon::Ribbon, ::Val{:Strand}, segment_range::UnitRange{Int}, chain::Protein.Chain)
    ca_points = Protein.alphacarbon_coords(chain)
    o_points = Protein.oxygen_coords(chain) .|> eltype(ca_points)
    return strand_surface(ribbon.attributes, ca_points, o_points; segment_range)
end

function render!(ribbon::Ribbon, chain::Protein.Chain, colors)
    for (ss_name, segment_range) in segments(chain)
        surface_vertices = get_surface_vertices(ribbon, Val(ss_name), segment_range, chain)
        surface!(ribbon, eachslice(surface_vertices, dims=1)...,
            color=expand_colors(colors[segment_range], size(surface_vertices, 2)), colormap=ribbon.colormap, colorrange=(0,1))
    end
    return ribbon
end

function draw_lines_between_subchains!(ribbon::Ribbon, subchains::AbstractVector{Protein.Chain}, color; linewidth=2)
    for (i, j) in zip(1:length(subchains)-1, 2:length(subchains))
        startpoint, endpoint = subchains[i].backbone[end], subchains[j].backbone[begin]
        n_segments = trunc(Int, norm(endpoint - startpoint) / 0.8)
        xs, ys, zs = [LinRange(startpoint[i], endpoint[i], 2*n_segments) for i in 1:3]
        linesegments!(ribbon, xs, ys, zs; linewidth=linewidth, color=color, transparency=true)
    end
    return ribbon
end

function render!(ribbon::Ribbon, chains::Vector{Protein.Chain})
    for (chain, colors) in zip(chains, ribbon.colors[])
        subchain_ranges = split(chain)
        if length(subchain_ranges) == 1
            render!(ribbon, chain, colors)
        else
            subchains = [chain[r] for r in subchain_ranges]
            for (subchain, r) in zip(subchains, subchain_ranges)
                render!(ribbon, subchain, colors[r])
            end
            draw_lines_between_subchains!(ribbon, subchains, :lightgray)
        end
    end
    return ribbon
end