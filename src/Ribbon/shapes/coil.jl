function coil_surface(
    attributes::Attributes,
    all_ca_points::AbstractMatrix{T};
    segment_range::UnitRange{Int},
) where T <: Real
    diameter = attributes.coil_diameter[]
    spline_quality = attributes.coil_spline_quality[]
    slice_quality = attributes.coil_slice_quality[]

    n_path_points = length(segment_range) * spline_quality
    path = spline(all_ca_points; N=n_path_points, r=segment_range)

    tangents = stack(path_tangent(path[:, max(1, i-1)], path[:, i], path[:, min(n_path_points, i+1)]) for i in 1:n_path_points)
    normals = zeros(T, 3, n_path_points)

    normals[:, 1] = abs(tangents[1, 1]) < 0.9 ? [1, 0, 0] : [0, 1, 0]
    for i in 2:n_path_points
        prev_normal = normals[:, i-1]
        t = tangents[:, i]
        projected_normal = prev_normal - dot(t, prev_normal) * t
        if norm(projected_normal) < 1e-5
            projected_normal = cross(tangents[:, i - 1], t)
        end
        normals[:, i] = normalize!(projected_normal)
    end

    surface_vertices = zeros(T, 3, n_path_points, slice_quality)
    for i in 1:n_path_points
        t = normalize(tangents[:, i])
        n = normalize(normals[:, i])
        b = cross(n, t)
        for (j, angle) in enumerate(range(0, 2π, slice_quality))
            offset = 0.5 * diameter * (cos(angle) * n + sin(angle) * b)
            surface_vertices[:, i, j] = path[:, i] + offset
        end
    end
    return surface_vertices
end

function get_surface_segment(ribbon::Ribbon, ::Val{:Coil}, segment_range::UnitRange{Int}, chain::Protein.Chain)
    adjusted_range = max(1, segment_range.start - 1):min(length(chain), segment_range.stop + 1)
    points = Protein.alphacarbon_coords(chain)
    return coil_surface(ribbon.attributes, points; segment_range=adjusted_range), adjusted_range
end