function coil_surface(
    attributes::Attributes,
    all_points::AbstractMatrix{T};
    segment_range::UnitRange{Int} = 1:size(all_midpoints, 2),
) where T <: Real
    radius = attributes.coil_radius[]
    spline_quality = attributes.coil_spline_quality[]
    slice_quality = attributes.coil_slice_quality[]

    points = all_points[:, segment_range]
    n_points = size(points, 2)
    n_path_points = n_points * spline_quality
    path = spline(all_points; N=n_path_points, r=segment_range)

    tangents = stack(path_tangent(path[:, max(1, i-1)], path[:, i], path[:, min(n_path_points, i+1)]) for i in 1:n_path_points)
    normals = zeros(T, 3, n_path_points)

    normals[:, 1] = abs(tangents[1, 1]) < 0.9 ? [1, 0, 0] : [0, 1, 0]
    for idx in 2:n_path_points
        prev_normal = normals[:, idx-1]
        t = tangents[:, idx]
        projected_normal = prev_normal - dot(t, prev_normal) * t
        if norm(projected_normal) < 1e-5
            projected_normal = cross(tangents[:, idx - 1], t)
        end
        normals[:, idx] = normalize!(projected_normal)
    end

    angles = LinRange(0, 2π, slice_quality)
    surface_vertices = zeros(T, 3, n_path_points, length(angles))
    for idx in 1:n_path_points
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