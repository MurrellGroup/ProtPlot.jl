# uses surrounding points to move the ends of the paths slightly away from each other, if they are equal
# such that the normals can be calculated correctly
function deintersect_ends!(path1::AbstractMatrix{T}, path2::AbstractMatrix{T}) where T <: Real
    if path1[:,1] == path2[:,1]
        next_binormal = path1[:,2] - path2[:,2]
        path1[:,1] += next_binormal .* 0.1
        path2[:,1] -= next_binormal .* 0.1
    end
    if path1[:,end] == path2[:,end]
        prev_binormal = path1[:,end-1] - path2[:,end-1]
        path1[:,end] += prev_binormal .* 0.1
        path2[:,end] -= prev_binormal .* 0.1
    end
end

arrow_function(l=0.5, w=0.5, W=1.0) = let hw = w/2, hW = W/2 
    t -> t > l ? hW*(t-1)/(l-1) : hw
end

function strand_midpoints(ca_points::AbstractMatrix{T}) where T <: Real
    N = size(ca_points, 2)
    midpoints = copy(ca_points)
    for idx in 2:N-1
        midpoints[:, idx] .= sum(ca_points[:, idx-1:idx+1], dims=2) / 3
    end
    return midpoints
end

function correct_alternating_normals!(normals::AbstractMatrix{T}) where T <: Real
    for (i, normal) in enumerate(eachcol(normals))
        if i % 2 == 0
            normal .= -normal
        end
    end
    return normals
end

function strand_surface(
    attributes::Attributes,
    all_ca_points::AbstractMatrix{T},
    all_o_points::AbstractMatrix{T};
    segment_range::UnitRange{Int} = 1:size(all_midpoints, 2),
) where T <: Real
    width = attributes.strand_width[]
    thickness = attributes.strand_thickness[]
    spline_quality = attributes.strand_spline_quality[]
    arrow_head_length = attributes.strand_arrow_head_length[]
    arrow_head_width = attributes.strand_arrow_head_width[]

    ca_points = all_ca_points[:, segment_range]
    o_points = all_o_points[:, segment_range]

    unnormalized_tangents = [ca_points[:, 1:end-1] - ca_points[:, 2:end];; ca_points[:, end-1] - ca_points[:, end]]
    binormaloids = o_points - ca_points

    normals = stack(normalize.(cross.(eachcol(unnormalized_tangents), eachcol(binormaloids))))
    correct_alternating_normals!(normals)

    binormals = stack(normalize.(cross.(eachcol(unnormalized_tangents), eachcol(normals))))

    midpoints = strand_midpoints(ca_points)

    n_points = size(midpoints, 2)
    n_path_points = n_points * spline_quality
    unnormalized_tangent_path = spline(unnormalized_tangents, N=n_path_points) ./ spline_quality
    normal_path = spline(normals, N=n_path_points)
    binormal_path = spline(binormals, N=n_path_points)
    midpoint_path = spline(midpoints, N=n_path_points)

    cumulative_path_length = cumsum(norm.(eachcol(unnormalized_tangent_path)))
    path_length = cumulative_path_length[end]
    arrow_body_length = path_length - arrow_head_length

    arrow_start_index = findfirst(>(arrow_body_length), cumulative_path_length)
    l = arrow_start_index / n_path_points
    arrow = arrow_function(l, width, arrow_head_width)

    surface_vertices = zeros(T, 3, n_path_points, 5)
    for i in 1:n_path_points
        t = (i - 1) / (n_path_points - 1)
        next_t = i / (n_path_points - 1)
        should_move_point = t < l < next_t # Ensure the arrow head has a 90-degree angle

        midpoint = midpoint_path[:, i + should_move_point]
        normal = normal_path[:, i + should_move_point]
        binormal = binormal_path[:, i + should_move_point]

        half_normal = normal * (thickness / 2)
        arrow_vector = binormal * arrow(t)
        surface_vertices[:, i, 1] = midpoint .+ half_normal .+ arrow_vector
        surface_vertices[:, i, 2] = midpoint .+ half_normal .- arrow_vector
        surface_vertices[:, i, 3] = midpoint .- half_normal .- arrow_vector
        surface_vertices[:, i, 4] = midpoint .- half_normal .+ arrow_vector
        surface_vertices[:, i, 5] = surface_vertices[:, i, 1]
    end

    surface_vertices[:, 1, :] .= midpoint_path[:, 1]

    return surface_vertices
end