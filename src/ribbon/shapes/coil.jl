function coil_surface(
    points::AbstractMatrix{T};
    radius = 1.0,
    spline_quality = 20,
    slice_quality = 20,
    ghost_control_start = nothing,
    ghost_control_end = nothing,
    kwargs...
) where T <: Real

    path = if !isnothing(ghost_control_start) && !isnothing(ghost_control_end)
        spline(points, ghost_control_start, ghost_control_end, m=spline_quality, k=min(3, size(points, 2)-1))
    else
        spline_quality = size(points, 2) == 2 ? 2 : spline_quality
        spline(points, m=spline_quality, k=min(3, size(points, 2)-1))
    end

    N = size(path, 2)
    angles = LinRange(0, 2π, slice_quality)
    surface_vertices = zeros(T, 3, N, length(angles))

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

        for (jdx, angle) in enumerate(angles)
            offset = radius * (cos(angle) * n + sin(angle) * b)
            surface_vertices[:, idx, jdx] = path[:, idx] + offset
        end
    end
    return surface_vertices
end