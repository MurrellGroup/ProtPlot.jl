include("shapes/shapes.jl")

function grid_to_mesh!(vertices::Vector{Point3f}, faces::Vector{GLTriangleFace}, vertex_colors::Vector{RGBAf},
                       surface_vertices::Array{<:Real,3}, path_colors::AbstractVector{RGBAf})
    M, N = size(surface_vertices, 2), size(surface_vertices, 3)
    offset = length(vertices)

    # Add vertices
    for i in 1:M, j in 1:N
        push!(vertices, Point3f(surface_vertices[1, i, j], surface_vertices[2, i, j], surface_vertices[3, i, j]))
    end

    # Add faces (triangulated quads)
    for i in 1:M-1, j in 1:N-1
        v1 = offset + (i-1)*N + j
        v2 = offset + i*N + j
        v3 = offset + i*N + j + 1
        v4 = offset + (i-1)*N + j + 1
        push!(faces, GLTriangleFace(v1, v2, v3))
        push!(faces, GLTriangleFace(v1, v3, v4))
    end

    # Add per-vertex colors (same color for all vertices at same path position)
    for i in 1:M
        c = path_colors[i]
        for j in 1:N
            push!(vertex_colors, c)
        end
    end

    return nothing
end

# Resolve scalar values through a colormap to produce RGBAf colors
function resolve_colormap(values::AbstractVector{<:Real}, colormap, colorrange)
    cmap = Makie.to_colormap(colormap)
    cmin, cmax = Float32.(colorrange)
    n = length(cmap)
    return map(values) do v
        t = clamp((Float32(v) - cmin) / (cmax - cmin), 0f0, 1f0)
        idx = 1f0 + t * (n - 1)
        i = clamp(floor(Int, idx), 1, n - 1)
        frac = idx - i
        c1 = RGBAf(cmap[i])
        c2 = RGBAf(cmap[min(i + 1, n)])
        RGBAf(
            c1.r * (1 - frac) + c2.r * frac,
            c1.g * (1 - frac) + c2.g * frac,
            c1.b * (1 - frac) + c2.b * frac,
            c1.alpha * (1 - frac) + c2.alpha * frac,
        )
    end
end

function add_segments!(vertices, faces, vertex_colors, ribbon, chain_backbone, secondary_structure, color)
    for (ss, segment_range) in segments(clean_secondary_structure!(secondary_structure))
        surface_vertices, adjusted_range = get_surface_segment(ribbon, Val(ss), segment_range, chain_backbone)
        expanded = vec(expand_colors(
            color[adjusted_range],
            size(surface_vertices, 2),
            [0.5; ones(length(adjusted_range)-2); 0.5]))
        grid_to_mesh!(vertices, faces, vertex_colors, surface_vertices, RGBAf.(expanded))
    end
end

function render_lines_between_subchains!(ribbon::Ribbon, chain_backbone_partitions::Vector{Array{T,3}}, color=:lightgray, linewidth=2) where T<:Real
    all_points = Point3f[]
    for i in 1:length(chain_backbone_partitions)-1
        startpoint, endpoint = @views chain_backbone_partitions[i][:, 2, end], chain_backbone_partitions[i+1][:, 2, 1]
        n_segments = max(1, trunc(Int, norm(endpoint - startpoint) / 0.8))
        coords = [range(startpoint[d], endpoint[d], 2*n_segments) for d in 1:3]
        for j in eachindex(coords[1])
            push!(all_points, Point3f(coords[1][j], coords[2][j], coords[3][j]))
        end
    end
    if !isempty(all_points)
        linesegments!(ribbon, all_points; linewidth=linewidth, color=color, transparency=true)
    end
    return nothing
end

function render!(ribbon::Ribbon, chain_backbones::Vector{Array{T,3}}) where T<:Real
    secondary_structures = ribbon.secondary_structures[]
    colors = ribbon.colors[]
    length(chain_backbones) == length(secondary_structures) || throw(ArgumentError("Chains and secondary structures vector must have the same length"))
    length(chain_backbones) == length(colors) || throw(ArgumentError("Chains and colors vector must have the same length"))

    # Determine color type and resolve all colors to RGBAf
    color_eltype = mapreduce(eltype, promote_type, colors)
    is_colorant = color_eltype <: Colorant

    vertices = Point3f[]
    faces = GLTriangleFace[]
    vertex_colors = RGBAf[]

    for (chain_backbone, secondary_structure, color) in zip(chain_backbones, secondary_structures, colors)
        chain_length = size(chain_backbone, 3)
        chain_length == length(secondary_structure) || throw(ArgumentError("coordinates and secondary structure size must match"))
        chain_length == length(color) || throw(ArgumentError("Chain and color must have the same length"))

        # Resolve to RGBAf: either from Colorants directly or through colormap
        resolved_color = if is_colorant
            RGBAf.(color)
        else
            resolve_colormap(color, ribbon.colormap[], ribbon.colorrange[])
        end

        partition_ranges = get_subchain_ranges(chain_backbone)
        chain_backbone_partitions = [chain_backbone[:, :, r] for r in partition_ranges]
        secondary_structure_partitions = [secondary_structure[r] for r in partition_ranges]
        color_partitions = [resolved_color[r] for r in partition_ranges]

        for (sub_backbone, sub_ss, sub_color) in zip(chain_backbone_partitions, secondary_structure_partitions, color_partitions)
            if size(sub_backbone, 3) > 1
                add_segments!(vertices, faces, vertex_colors, ribbon, sub_backbone, sub_ss, sub_color)
            else
                sphere = Sphere(Point3f(sub_backbone[:, 2, 1]), ribbon.coil_diameter[])
                mesh!(ribbon, sphere; color=:lightgray)
            end
        end
        ribbon.show_gaps[] && render_lines_between_subchains!(ribbon, chain_backbone_partitions)
    end

    if !isempty(vertices)
        gb_mesh = GeometryBasics.normal_mesh(
            GeometryBasics.Mesh(vertices, faces)
        )
        mesh!(ribbon, gb_mesh; color=vertex_colors)
    end

    return nothing
end
