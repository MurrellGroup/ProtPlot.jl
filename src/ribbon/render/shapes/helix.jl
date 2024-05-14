# normal depends on second derivative of the path
# does weird twist when second derivative changes too quickly
function helix_surface(points::AbstractMatrix{T}, attributes) where T <: Real
    radius = attributes.helix_radius
    width_factor = attributes.helix_width
    thickness_factor = attributes.helix_thickness
    spline_quality = attributes.helix_spline_quality
    slice_quality = attributes.helix_slice_quality

    path = spline(points, m=spline_quality, k=min(3, size(points, 2)-1))
    N = size(path, 2)
    angles = LinRange(0, 2π, slice_quality)
    surface_vertices = zeros(T, 3, N, length(angles))
    for idx in 1:N
        three_path_points = eachcol(@view(path[:, min(max(begin,idx-1),end-2):min(max(begin,idx-1)+2,end)])) # understandable expressions are for the weak
        t = path_tangent(three_path_points...)
        n = curved_path_normal(three_path_points...)
        n = cross(t, cross(n, t))

        b = cross(n, t)
        for (jdx, v) in enumerate(angles)
            offset = radius .* (thickness_factor*cos(v) .* n .+ width_factor*sin(v) .* b)
            surface_vertices[:, idx, jdx] = path[:, idx] .- offset
        end
    end
    return surface_vertices
end