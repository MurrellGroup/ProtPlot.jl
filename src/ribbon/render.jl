function render!(
    container,
    segment::Segment{Loop},
    colors::AbstractVector{<:RGB} = [colorant"red", colorant"yellow"];
    radius = 0.25,
    spline_quality = 20,
    slice_quality = 20,
    kwargs...
)
    extended_segment = extend_segment(segment, 0:length(segment)+1)
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    controls = @view alphacarbon_coord_matrix(segment.backbone)[:, (isone(end) ? 1 : 2):end]
    coords = [startpoint controls endpoint]
    surface_vertices = tube_surface(coords,
        radius=radius, spline_quality=spline_quality, slice_quality=slice_quality,
        ghost_control_start=segment_startpoint(extended_segment), ghost_control_end=segment_endpoint(extended_segment)
    )
    N = size(surface_vertices, 2)
    color_matrix = expand_colors(colors, N)
    surface!(container, eachslice(surface_vertices, dims=1)..., color=color_matrix)

    return container
end

function render!(
    container,
    segment::Segment{Helix},
    colors::AbstractVector{<:RGB} = [colorant"lime", colorant"cyan"];
    helix_radius = 1.0,
    helix_width = 1.0,
    helix_thickness = 0.25,
    spline_quality = 20,
    slice_quality = 20,
    kwargs...
)
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    controls = @view alphacarbon_coord_matrix(segment.backbone)[:, 2:end] # startpoint is first point instead. including first N *and* CA could mess with normals
    coords = hcat(startpoint, controls, endpoint)
    surface_vertices = helix_surface(coords,
        radius=helix_radius, width_factor=helix_width, thickness_factor=helix_thickness,
        spline_quality=spline_quality, slice_quality=slice_quality,
    )
    N = size(surface_vertices, 2)
    color_matrix = expand_colors(colors, N)
    surface!(container, eachslice(surface_vertices, dims=1)..., color=color_matrix)

    return container
end

function render!(
    container,
    segment::Segment{Strand},
    colors::AbstractVector{<:RGB} = [colorant"blue", colorant"magenta"];
    strand_width = 2.0,
    strand_thickness = 0.5,
    spline_quality = 20,
    kwargs...
)
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    oxygen_coords_side1 = @view oxygen_coord_matrix(segment.backbone)[:, 1:2:end-1]
    oxygen_coords_side2 = @view oxygen_coord_matrix(segment.backbone)[:, 2:2:end]
    coords1 = hcat(startpoint, oxygen_coords_side1, endpoint)
    coords2 = hcat(startpoint, oxygen_coords_side2, endpoint)
    surface_vertices = arrow_surface(coords1, coords2,
        width=strand_width, thickness=strand_thickness, spline_quality=spline_quality,
    )
    N = size(surface_vertices, 2)
    color_matrix = expand_colors(colors, N) 
    surface!(container, eachslice(surface_vertices, dims=1)..., color=color_matrix)

    return container
end

function render!(
    container,
    chain::Chain,
    colors::AbstractVector{<:RGB};
    kwargs...
)
    @assert length(chain) == length(colors)
    @assert !has_missing_ss(chain)
    for segment in segments(chain)
        render!(container, segment, colors[segment.range]; kwargs...)
    end

    return container
end

function render!(
    container,
    protein::Protein;
    colorscheme::Union{ColorScheme, Symbol} = :jet,
    color_vectors::AbstractVector{C} = [LinRange(0, 1, length(chain)) for chain in protein],
    kwargs...
) where C <: AbstractVector{<:Union{Real, RGB}}
    colorscheme isa Symbol && (colorscheme = colorschemes[colorscheme])
    if eltype(C) <: Real
        color_vectors = [colorscheme[color_vector] for color_vector in color_vectors]
    end
    has_missing_ss(protein) && assign_secondary_structure!(protein)
    remove_singleton_strands!.(protein) # TODO: don't mutate
    @assert length(protein) == length(color_vectors)
    @assert length.(protein) == length.(color_vectors)
    for (chain, colors) in zip(protein, color_vectors)
        render!(container, chain, colors; kwargs...)
    end

    return container
end