path_tangent(p1, p2, p3) = normalize!(p3 - p1)
curved_path_normal(p1, p2, p3) = -normalize!(p1 - p2 + p3 - p2)

function tube_surface(
    points::AbstractMatrix{T},
    radius::Real = 0.5;
    spline_quality = 10,
    tube_quality = 20,
    ghost_control_start = nothing,
    ghost_control_end = nothing,
) where T <: Real

    path = if !isnothing(ghost_control_start) && !isnothing(ghost_control_end)
        spline(points, ghost_control_start, ghost_control_end, m=spline_quality, k=min(3, size(points, 2)-1))
    else
    spline_quality = size(points, 2) == 2 ? 2 : spline_quality
        spline(points, m=spline_quality, k=min(3, size(points, 2)-1))
    end

    N = size(path, 2)
    tube_angles = LinRange(0, 2π, tube_quality)
    surface_vertices = zeros(T, 3, N, length(tube_angles))

    # Precompute tangents and initialize normals
    tangents = stack(path_tangent(path[:, max(1, i-1)], path[:, i], path[:, min(N, i+1)]) for i in 1:N)
    normals = zeros(T, 3, N)

    # Initial normal vector (arbitrary, but not aligned with the first tangent)
    normals[:, 1] = abs(tangents[1, 1]) < 0.9 ? [1, 0, 0] : [0, 1, 0]

    # Propagate the normal vector along the path
    for idx in 2:N
        prev_normal = normals[:, idx-1]
        t = tangents[:, idx]
        projected_normal = prev_normal - dot(t, prev_normal) * t
        if norm(projected_normal) < 1e-5
            projected_normal = cross(tangents[:, idx - 1], t)
        end
        normals[:, idx] = normalize!(projected_normal)
    end

    # Generate surface vertices
    for idx in 1:N
        t = tangents[:, idx]
        n = normals[:, idx]
        b = cross(n, t)

        for (jdx, angle) in enumerate(tube_angles)
            offset = radius * (cos(angle) * n + sin(angle) * b)
            surface_vertices[:, idx, jdx] = path[:, idx] + offset
        end
    end
    return surface_vertices
end

# Function to generate a tube surface given a path and radius
# normal depends on second derivative of the path
# does weird twist when second derivative changes too quickly
function helix_surface(
    points::AbstractMatrix{T},
    radius::Real = 0.5;
    spline_quality = 10,
    tube_quality = 20,
    x_elongation = 1,
    z_elongation = 1,
) where T <: Real
    path = spline(points, m=spline_quality, k=min(3, size(points, 2)-1))
    N = size(path, 2)
    tube_angles = LinRange(0, 2π, tube_quality)
    surface_vertices = zeros(T, 3, N, length(tube_angles))
    for idx in 1:N
        three_path_points = eachcol(@view(path[:, min(max(begin,idx-1),end-2):min(max(begin,idx-1)+2,end)])) # understandable expressions are for the weak
        t = path_tangent(three_path_points...)
        n = curved_path_normal(three_path_points...)
        n = cross(t, cross(n, t))

        b = cross(n, t)
        for (jdx, v) in enumerate(tube_angles)
            offset = radius .* (x_elongation*cos(v) .* n .+ z_elongation*sin(v) .* b)
            surface_vertices[:, idx, jdx] = path[:, idx] .- offset
        end
    end
    return surface_vertices
end

# uses surrounding points to move the ends of the paths slightly away from each other, if they are equal
# such that the normals can be calculated correctly
function deintersect_ends!(path1::AbstractMatrix{T}, path2::AbstractMatrix{T}) where T <: Real
    @assert size(path1, 2) == size(path2, 2)
    @assert size(path1, 1) >= 3
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

# 0 <= t <= 1
# l is the length of the arrow body
# w is the width of the arrow body
# W is the width of the arrow head
# the length of the arrow head becomes 1 - l
# the shape is normalized such that length of the entire arrow is 1
arrow_function(l=0.5, w=0.5, W=1.0) = t -> t > l ? W*(t-1)/(l-1) : w

# TODO: just use oxygen atoms for calculating paths

function arrow_surface(
    points1::AbstractMatrix{T},
    points2::AbstractMatrix{T};
    width = 1.0,
    thickness = 0.3,
    spline_quality = 10,
) where T <: Real
    half_thickness = thickness / 2
    max_L = max(size(points1, 2), size(points2, 2))
    path1 = spline(points1, N=max_L*spline_quality, k=min(3, size(points1, 2)-1))
    path2 = spline(points2, N=max_L*spline_quality, k=min(3, size(points2, 2)-1))
    deintersect_ends!(path1, path2) # this is a hack to make the normals work
    @assert size(path1, 2) == size(path2, 2)
    N = size(path1, 2)

    midpath = (path1 + path2) / 2
    unnormalized_tangents = hcat(midpath[:, 1:end-1] - midpath[:, 2:end], midpath[:, end-1] - midpath[:, end])
    tangents = mapslices(normalize, unnormalized_tangents, dims=1)
    almost_binormals = mapslices(normalize, midpath - path1, dims=1)
    normals = stack(normalize!.(cross.(eachcol(tangents), eachcol(almost_binormals))))
    binormals = stack(cross.(eachcol(normals), eachcol(tangents)))

    cumulative_length_of_path = cumsum(norm.(eachcol(unnormalized_tangents)))
    length_of_path = cumulative_length_of_path[end]
    arrow_head_length = 3.0
    arrow_body_length = length_of_path - arrow_head_length
    l = findfirst(>(arrow_body_length), cumulative_length_of_path) / N
    arrow = arrow_function(l, width, width*2.0)

    surface_vertices = zeros(T, 3, N, 5)
    for idx in 1:N
        midpoint = midpath[:, idx]
        normal = normals[:, idx]
        binormal = binormals[:, idx]

        half_normal = half_thickness .* normal
        arrow_vector = binormal * arrow((idx-1)/(N-1))
        surface_vertices[:, idx, 1] = midpoint .+ half_normal .+ arrow_vector
        surface_vertices[:, idx, 2] = midpoint .+ half_normal .- arrow_vector
        surface_vertices[:, idx, 3] = midpoint .- half_normal .- arrow_vector
        surface_vertices[:, idx, 4] = midpoint .- half_normal .+ arrow_vector
        surface_vertices[:, idx, 5] = surface_vertices[:, idx, 1]
    end

    surface_vertices[:, 1, :] .= midpath[:, 1]

    return surface_vertices
end