path_tangent(p1, p2, p3) = normalize!(p3 - p1)
curved_path_normal(p1, p2, p3) = -normalize!(p1 - p2 + p3 - p2)

nitrogen_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 1)
alphacarbon_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 2)
carbon_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 3)
oxygen_coord_matrix(chain::AbstractChain) = atom_coord_matrix(chain, 4)

const DEFAULT_COLORSCHEME = :jet

function tube_surface(
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
    tube_surface_vertices = zeros(T, 3, N, length(tube_angles))
    for idx in 1:N
        three_path_points = eachcol(@view(path[:, min(max(begin,idx-1),end-2):min(max(begin,idx-1)+2,end)])) # understandable expressions are for the weak
        t = path_tangent(three_path_points...)
        n = curved_path_normal(three_path_points...)
        n = cross(t, cross(n, t))

        b = cross(n, t)        
        for (jdx, v) in enumerate(tube_angles)
            offset = radius .* (x_elongation*cos(v) .* n .+ z_elongation*sin(v) .* b)
            tube_surface_vertices[:, idx, jdx] = path[:, idx] .- offset
        end
    end
    return tube_surface_vertices
end

# Function to generate a tube surface given a path and radius
# normal depends on second derivative of the path
# does weird twist when second derivative changes too quickly
function tube(
    points::AbstractMatrix{T},
    radius::Real = 0.5;
    spline_quality = 10,
    tube_quality = 20,
    x_elongation = 1,
    z_elongation = 1,
    color_start = 0,
    color_end = 1,
    colorscheme::Union{ColorScheme, Symbol} = DEFAULT_COLORSCHEME,
) where T <: Real
    surface_vertices = tube_surface(
        points, radius,
        spline_quality=spline_quality, tube_quality=tube_quality,
        x_elongation=x_elongation, z_elongation=z_elongation)
    colorscheme = colorscheme isa ColorScheme ? colorscheme : colorschemes[colorscheme]
    color_matrix = repeat(colorscheme[LinRange(color_start, color_end, size(surface_vertices, 2))], inner=(spline_quality, 1))
    return surface_vertices, color_matrix
end

# uses surrounding points to move the ends of the paths slightly away from each other, if they are equal
# such that the normals can be calculated correctly
function deintersect_ends!(path1::AbstractMatrix{T}, path2::AbstractMatrix{T}) where T <: Real
    @assert size(path1, 2) == size(path2, 2)
    @assert size(path1, 1) >= 3
    if path1[:,1] == path2[:,1]
        next_binormal = path1[:,2] - path2[:,2]
        path1[:,1] += next_binormal .* 0.01
        path2[:,1] -= next_binormal .* 0.01
    end
    if path1[:,end] == path2[:,end]
        prev_binormal = path1[:,end-1] - path2[:,end-1]
        path1[:,end] += prev_binormal .* 0.01
        path2[:,end] -= prev_binormal .* 0.01
    end
end

function sheet_surface(
    points1::AbstractMatrix{T},
    points2::AbstractMatrix{T};
    thickness = 0.5,
    width_pad = 0.0,
    spline_quality = 4,
) where T <: Real
    half_thickness = thickness / 2
    max_L = max(size(points1, 2), size(points2, 2))
    path1 = spline(points1, N=max_L*spline_quality, k=min(3, size(points1, 2)-1))
    path2 = spline(points2, N=max_L*spline_quality, k=min(3, size(points2, 2)-1))
    deintersect_ends!(path1, path2)
    @assert size(path1, 2) == size(path2, 2)
    N = size(path1, 2)
    sheet_surface_vertices = zeros(T, 3, N, 5)
    for idx in 1:N
        point1 = path1[:, idx]
        point2 = path2[:, idx]
        three_path1_points = eachcol(@view(path1[:, min(max(begin,idx-1),end-2):min(max(begin,idx-1)+2,end)]))
        three_path2_points = eachcol(@view(path2[:, min(max(begin,idx-1),end-2):min(max(begin,idx-1)+2,end)]))

        tangent1 = path_tangent(three_path1_points...)
        binormal1 = normalize!(point1 - point2)
        normal1 = normalize!(cross(binormal1, tangent1))

        tangent2 = path_tangent(three_path2_points...)
        binormal2 = normalize!(point2 - point1)
        normal2 = normalize!(cross(binormal2, tangent2))

        padding = width_pad .* binormal1

        sheet_surface_vertices[:, idx, 1] = point1 .+ half_thickness .* normal1 .+ padding
        sheet_surface_vertices[:, idx, 2] = point2 .- half_thickness .* normal2 .- padding
        sheet_surface_vertices[:, idx, 3] = point2 .+ half_thickness .* normal2 .- padding
        sheet_surface_vertices[:, idx, 4] = point1 .- half_thickness .* normal1 .+ padding
        sheet_surface_vertices[:, idx, 5] = sheet_surface_vertices[:, idx, 1]
    end
    return sheet_surface_vertices
end

function sheet(
    points1::AbstractMatrix{T},
    points2::AbstractMatrix{T};
    thickness = 0.5,
    width_pad = 0.0,
    spline_quality = 4,
    color_start = 0,
    color_end = 1,
    colorscheme::Union{ColorScheme, Symbol} = DEFAULT_COLORSCHEME,
) where T <: Real
    surface_vertices = sheet_surface(
        points1, points2,
        thickness=thickness, width_pad=width_pad, spline_quality=spline_quality)
    colorscheme = colorscheme isa ColorScheme ? colorscheme : colorschemes[colorscheme]
    color_matrix = repeat(colorscheme[LinRange(color_start, color_end, size(surface_vertices, 2))], inner=(spline_quality, 1))
    return surface_vertices, color_matrix
end