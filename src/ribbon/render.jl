export ribbon, ribbon!

function render!(
    container,
    segment::Segment{Loop},
    colors::AbstractVector{<:RGB} = [colorant"red", colorant"yellow"],
)
    extended_segment = extend_segment(segment, 0:length(segment)+1)
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    controls = @view alphacarbon_coord_matrix(segment.backbone)[:, (isone(end) ? 1 : 2):end]
    coords = [startpoint controls endpoint]
    surface_vertices = tube_surface(coords, 0.25, spline_quality=20, tube_quality=20,
        ghost_control_start=segment_startpoint(extended_segment), ghost_control_end=segment_endpoint(extended_segment))
    N = size(surface_vertices, 2)
    color_matrix = expand_colors(colors, N)
    surface!(container, eachslice(surface_vertices, dims=1)..., color=color_matrix)
end

function render!(
    container,
    segment::Segment{Helix},
    colors::AbstractVector{<:RGB} = [colorant"lime", colorant"cyan"],
)
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    controls = @view alphacarbon_coord_matrix(segment.backbone)[:, 2:end] # startpoint is first point instead. including first N *and* CA could mess with normals
    coords = hcat(startpoint, controls, endpoint)
    surface_vertices = helix_surface(coords, 1, spline_quality=20, tube_quality=20, x_elongation=0.25)
    N = size(surface_vertices, 2)
    color_matrix = expand_colors(colors, N)
    surface!(container, eachslice(surface_vertices, dims=1)..., color=color_matrix)
end

function render!(
    container,
    segment::Segment{Strand},
    colors::AbstractVector{<:RGB} = [colorant"blue", colorant"magenta"],
)
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    oxygen_coords_side1 = @view oxygen_coord_matrix(segment.backbone)[:, 1:2:end-1]
    oxygen_coords_side2 = @view oxygen_coord_matrix(segment.backbone)[:, 2:2:end]
    coords1 = hcat(startpoint, oxygen_coords_side1, endpoint)
    coords2 = hcat(startpoint, oxygen_coords_side2, endpoint)
    surface_vertices = arrow_surface(coords1, coords2, width=1.1, thickness=0.5, spline_quality=20)
    N = size(surface_vertices, 2)
    color_matrix = expand_colors(colors, N) 
    surface!(container, eachslice(surface_vertices, dims=1)..., color=color_matrix)
end

function render!(
    container,
    chain::Chain,
    colors::AbstractVector{<:RGB},
)
    @assert length(chain) == length(colors)
    @assert !has_missing_ss(chain)
    for segment in segments(chain)
        render!(container, segment, colors[segment.range])
    end
end

function ribbon!(
    container,
    protein::Protein,
    color_vectors::AbstractVector{<:AbstractVector{<:RGB}},
)
    has_missing_ss(protein) && assign_secondary_structure!(protein)
    remove_singleton_strands!.(protein) # TODO: don't mutate
    @assert length(protein) == length(color_vectors)
    @assert length.(protein) == length.(color_vectors)
    for (chain, colors) in zip(protein, color_vectors)
        render!(container, chain, colors)
    end
end

function ribbon!(
    container,
    protein::Protein,
    float_color_vectors::AbstractVector{<:AbstractVector{<:AbstractFloat}},
    colorscheme::ColorScheme = ColorSchemes.jet,
)
    color_vectors = [colorscheme[float_color_vector] for float_color_vector in float_color_vectors]
    ribbon!(container, protein, color_vectors)
end

function ribbon!(
    container,
    protein::Protein,
    colorscheme::ColorScheme = ColorSchemes.jet,
)
    color_vectors = [smooth_color_vector(colorscheme, length(chain)) for chain in protein]
    ribbon!(container, protein, color_vectors)
end

function ribbon(
    protein::Protein,
    color_vectors::ColorScheme = [smooth_color_vector(colorscheme, length(chain)) for chain in protein],
)
    scene = Scene(backgroundcolor=:black)
    cam3d!(scene)
    ribbon!(scene, protein, color_vectors)
    center!(scene)
    display(scene)
end
