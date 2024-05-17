curved_path_normal(p1, p2, p3) = -normalize!(p1 - p2 + p3 - p2)

# normal depends on second derivative of the path
# does weird twist when second derivative changes too quickly
function helix_surface(
    attributes::Attributes,
    all_ca_points::AbstractMatrix{T};
    segment_range::UnitRange{Int} = 1:size(all_midpoints, 2)
) where T <: Real
    width = attributes.helix_width[]
    thickness = attributes.helix_thickness[]
    spline_quality = attributes.helix_spline_quality[]
    slice_quality = attributes.helix_slice_quality[]

    n_path_points = length(segment_range) * spline_quality
    path = spline(all_ca_points; N=n_path_points, r=segment_range)

    surface_vertices = zeros(T, 3, n_path_points, slice_quality)
    for i in 1:n_path_points
        three_path_points = eachcol(@view(path[:, min(max(begin,i-1),end-2):min(max(begin,i-1)+2,end)])) # understandable expressions are for the weak
        t = normalize(path_tangent(three_path_points...))
        n = normalize(cross(t, cross(curved_path_normal(three_path_points...), t)))
        b = cross(n, t)
        for (j, v) in enumerate(LinRange(0, 2π, slice_quality))
            offset = 0.5*thickness*cos(v) .* n .+ 0.5*width*sin(v) .* b
            surface_vertices[:, i, j] = path[:, i] .- offset
        end
    end
    
    surface_vertices[:, 1, :] .= path[:, 1]
    surface_vertices[:, end, :] .= path[:, end]

    return surface_vertices
end