include("shapes/shapes.jl")

function render!(ribbon::Ribbon, chain_backbone::Array{T,3}, secondary_structure::Vector{Int}, color::AbstractVector{<:Real}, colormap) where T<:Real
    for (ss, segment_range) in segments(clean_secondary_structure!(secondary_structure))
        surface_vertices, adjusted_range = get_surface_segment(ribbon, Val(ss), segment_range, chain_backbone)
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

function render!(ribbon::Ribbon, chain_backbone::Array{T,3}, secondary_structure::Vector{Int}, color::AbstractVector{<:Colorant}, colormap) where T<:Real
    render!(ribbon, chain_backbone, secondary_structure, range(0, 1, length(color)), color)
end

function render_lines_between_subchains!(ribbon::Ribbon, chain_backbone_partitions::Vector{Array{T,3}}, color=:lightgray, linewidth=2) where T<:Real
    for i in 1:length(chain_backbone_partitions)-1
        startpoint, endpoint = @views chain_backbone_partitions[i][:, 2, end], chain_backbone_partitions[i+1][:, 2, 1]
        n_segments = trunc(Int, norm(endpoint - startpoint) / 0.8)
        xs, ys, zs = [range(startpoint[i], endpoint[i], 2*n_segments) for i in 1:3]
        linesegments!(ribbon, xs, ys, zs; linewidth=linewidth, color=color, transparency=true)
    end
    return nothing
end

function render!(ribbon::Ribbon, chain_backbones::Vector{Array{T,3}}) where T<:Real
    secondary_structures = ribbon.secondary_structures[]
    colors = ribbon.colors[]
    length(chain_backbones) == length(secondary_structures) || throw(ArgumentError("Chains and secondary structures vector must have the same length"))
    length(chain_backbones) == length(colors) || throw(ArgumentError("Chains and colors vector must have the same length"))
    for (chain_backbone, secondary_structure, color) in zip(chain_backbones, secondary_structures, colors)
        chain_length = size(chain_backbone, 3)
        chain_length == length(secondary_structure) || throw(ArgumentError("coordinates and secondary structure size must match"))
        chain_length == length(color) || throw(ArgumentError("Chain and color must have the same length"))

        partition_ranges = get_subchain_ranges(chain_backbone)
        chain_backbone_partitions = [chain_backbone[:, :, r] for r in partition_ranges]
        secondary_structure_partitions = [secondary_structure[r] for r in partition_ranges]
        color_partitions = [color[r] for r in partition_ranges]
        for (chain_backbone, secondary_structure, color) in zip(chain_backbone_partitions, secondary_structure_partitions, color_partitions)
            if size(chain_backbone, 3) > 1
                render!(ribbon, chain_backbone, secondary_structure, color, ribbon.colormap)
            else
                sphere = Sphere(Point3f(chain_backbone[:, 2, 1]), ribbon.coil_diameter[])
                mesh!(ribbon, sphere; color=:lightgray)
            end
        end
        render_lines_between_subchains!(ribbon, chain_backbone_partitions)
    end
    return nothing
end
