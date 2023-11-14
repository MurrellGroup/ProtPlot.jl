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