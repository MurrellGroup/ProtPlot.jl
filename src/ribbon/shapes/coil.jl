function coil_surface(
    attributes::Attributes,
    all_points::AbstractMatrix{T};
    segment_range::UnitRange{Int} = 1:size(all_midpoints, 2),
) where T <: Real
    radius = attributes.coil_radius[]
    spline_quality = attributes.coil_spline_quality[]
    slice_quality = attributes.coil_slice_quality[]

    points = all_points[:, segment_range]

    spline_quality = size(points, 2) == 2 ? 2 : spline_quality
    path = spline(all_points; N=size(points, 2) * spline_quality, r=segment_range)

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