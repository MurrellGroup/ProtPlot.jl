# uses surrounding points to move the ends of the paths slightly away from each other, if they are equal
# such that the normals can be calculated correctly
function deintersect_ends!(path1::AbstractMatrix{T}, path2::AbstractMatrix{T}) where T <: Real
    size(path1, 2) == size(path2, 2) || throw(ArgumentError("The paths must have the same number of points."))
    size(path1, 1) >= 3 || throw(ArgumentError("The paths must have at least 3 dimensions."))
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
arrow_function(l=0.5, w=0.5, W=1.0) = let hw = w/2, hW = W/2 
    t -> t > l ? hW*(t-1)/(l-1) : hw
end

function arrow_surface(points1::AbstractMatrix{T}, points2::AbstractMatrix{T}, attributes) where T <: Real
    width = attributes.strand_width
    thickness = attributes.strand_thickness
    arrow_head_length = attributes.strand_arrow_head_length
    arrow_head_width = attributes.strand_arrow_head_width
    spline_quality = attributes.strand_spline_quality

    if size(points1, 2) > 3 && size(points2, 2) > 3
        points1 = points1[:, [1:end-2; end]]
        points2 = points2[:, [1:end-2; end]]
    end
    half_thickness = thickness / 2
    max_L = max(size(points1, 2), size(points2, 2))
    path1 = spline(points1, N=max_L*spline_quality, k=min(3, size(points1, 2)-1))
    path2 = spline(points2, N=max_L*spline_quality, k=min(3, size(points2, 2)-1))
    deintersect_ends!(path1, path2) # this is a hack to make the normals work

    N = size(path1, 2)

    midpath = (path1 + path2) / 2
    unnormalized_tangents = hcat(midpath[:, 1:end-1] - midpath[:, 2:end], midpath[:, end-1] - midpath[:, end])
    tangents = mapslices(normalize, unnormalized_tangents, dims=1)
    almost_binormals = mapslices(normalize, midpath - path1, dims=1)
    normals = stack(normalize!.(cross.(eachcol(tangents), eachcol(almost_binormals))))
    binormals = stack(cross.(eachcol(normals), eachcol(tangents)))

    cumulative_length_of_path = cumsum(norm.(eachcol(unnormalized_tangents)))
    length_of_path = cumulative_length_of_path[end]
    arrow_body_length = length_of_path - arrow_head_length
    l = findfirst(>(arrow_body_length), cumulative_length_of_path) / N
    arrow = arrow_function(l, width, arrow_head_width)

    surface_vertices = zeros(T, 3, N, 5)
    for idx in 1:N
        t = (idx-1)/(N-1)
        next_t = idx/(N-1)
        should_move_point = t < l < next_t # hack to get the arrow head to have a 90 degree angle

        midpoint = midpath[:, idx + should_move_point]
        normal = normals[:, idx + should_move_point]
        binormal = binormals[:, idx + should_move_point]

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