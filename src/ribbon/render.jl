include("shapes/shapes.jl")

function render!(ribbon::Ribbon, chain::Protein.Chain, color::AbstractVector{<:Real}, colormap)
    for (ss_name, segment_range) in segments(chain)
        surface_vertices, adjusted_range = get_surface_segment(ribbon, Val(ss_name), segment_range, chain)
        surface!(ribbon, eachslice(surface_vertices, dims=1)...,
            color = expand_colors(
                color[adjusted_range],
                size(surface_vertices, 2),
                [0.5; ones(length(adjusted_range)-2); 0.5]),
            colormap = colormap,
            colorrange = (0, 1))
    end
    return nothing
end

function render!(ribbon::Ribbon, chain::Protein.Chain, color::AbstractVector{<:Colorant}, colormap)
    render!(ribbon, chain, range(0, 1, length(color)), color)
end

function render_lines_between_subchains!(ribbon::Ribbon, subchains::Vector{Protein.Chain}, color=:lightgray, linewidth=2)
    for i in 1:length(subchains)-1
        startpoint, endpoint = subchains[i].backbone[end-1], subchains[i+1].backbone[begin+1]
        n_segments = trunc(Int, norm(endpoint - startpoint) / 0.8)
        xs, ys, zs = [range(startpoint[i], endpoint[i], 2*n_segments) for i in 1:3]
        linesegments!(ribbon, xs, ys, zs; linewidth=linewidth, color=color, transparency=true)
    end
    return nothing
end

function render!(ribbon::Ribbon, chains::Vector{Protein.Chain})
    colors = ribbon.colors[]
    length(chains) == length(colors) || throw(ArgumentError("Chains and colors must have the same length"))
    for (chain, color) in zip(chains, colors)
        length(chain) > 1 || throw(ArgumentError("Chain must have at least 2 residues"))
        length(chain) == length(color) || throw(ArgumentError("Chain and color must have the same length"))
        subchain_ranges = get_subchain_ranges(chain)
        subchains = [chain[r] for r in subchain_ranges]
        subcolors = [color[r] for r in subchain_ranges]
        for (subchain, subcolor) in zip(subchains, subcolors)
            render!(ribbon, subchain, subcolor, ribbon.colormap)
        end
        render_lines_between_subchains!(ribbon, subchains)
    end
    return nothing
end