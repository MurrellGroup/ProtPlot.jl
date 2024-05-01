function render!(
    container,
    segment::Segment{Coil},
    colors::AbstractVector{<:RGB} = [colorant"red", colorant"yellow"];
    radius = 0.2,
    spline_quality = 20,
    slice_quality = 20,
    plots = nothing,
    kwargs...
)
    controls = Protein.alphacarbon_coords(Protein.Chain(segment.backbone))[:, (isone(end) ? 1 : 2):(end > 2 ? end-1 : end)]
    coords = [segment_startpoint(segment) controls segment_endpoint(segment)]
    surface_vertices = coil_surface(coords;
        radius=radius, spline_quality=spline_quality, slice_quality=slice_quality,
        ghost_control_start = segment.range.start == 1 ? nothing : segment.chain.backbone[3*segment.range.start-2],
        ghost_control_end = segment.range.stop == length(segment.chain) ? nothing : segment.chain.backbone[3*segment.range.stop+2],
        kwargs...
    )
    N = size(surface_vertices, 2)
    color_matrix = expand_colors(colors, N)
    p = surface!(container, eachslice(surface_vertices, dims=1)..., color=color_matrix)
    !isnothing(plots) && push!(plots, p)

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
    plots = nothing,
    kwargs...
)
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    controls = Protein.alphacarbon_coords(Protein.Chain(segment.backbone))[:, 2:end] # startpoint is first point instead. including first N *and* CA could mess with normals
    coords = hcat(startpoint, controls, endpoint)
    surface_vertices = helix_surface(coords;
        radius=helix_radius, width_factor=helix_width, thickness_factor=helix_thickness,
        spline_quality=spline_quality, slice_quality=slice_quality,
        kwargs...
    )
    N = size(surface_vertices, 2)
    color_matrix = expand_colors(colors, N)
    p = surface!(container, eachslice(surface_vertices, dims=1)..., color=color_matrix)
    !isnothing(plots) && push!(plots, p)

    return container
end

function render!(
    container,
    segment::Segment{Strand},
    colors::AbstractVector{<:RGB} = [colorant"blue", colorant"magenta"];
    strand_width = 2.0,
    strand_thickness = 0.5,
    spline_quality = 20,
    plots = nothing,
    kwargs...
)
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    oxygen_coords_side1 = @view Protein.oxygen_coords(Protein.Chain(segment.backbone))[:, 1:2:end-1]
    oxygen_coords_side2 = @view Protein.oxygen_coords(Protein.Chain(segment.backbone))[:, 2:2:end]
    coords1 = hcat(startpoint, oxygen_coords_side1, endpoint)
    coords2 = hcat(startpoint, oxygen_coords_side2, endpoint)
    surface_vertices = arrow_surface(coords1, coords2;
        width=strand_width, thickness=strand_thickness, spline_quality=spline_quality,
        kwargs...
    )
    N = size(surface_vertices, 2)
    color_matrix = expand_colors(colors, N) 
    p = surface!(container, eachslice(surface_vertices, dims=1)..., color=color_matrix)
    !isnothing(plots) && push!(plots, p)

    return container
end

function render!(
    container,
    chain::Protein.Chain,
    colors::AbstractVector{<:RGB};
    kwargs...
)
    @assert length(chain) == length(colors)
    @assert Protein.has_assigned_ss(chain)
    for segment in segments(chain)
        render!(container, segment, colors[segment.range]; kwargs...)
    end

    return container
end

function draw_lines_between_subchains!(container, subchains::AbstractVector{Protein.Chain}, color::RGB; linewidth=2, plots=nothing, kwargs...)
    for (i, j) in zip(1:length(subchains)-1, 2:length(subchains))
        startpoint, endpoint = subchains[i].backbone[end], subchains[j].backbone[begin]
        n_segments = trunc(Int, norm(endpoint - startpoint) / 0.8)
        xs, ys, zs = [LinRange(startpoint[i], endpoint[i], 2*n_segments) for i in 1:3]
        p = linesegments!(container, xs, ys, zs; linewidth=linewidth, color=color, transparency=true)
        !isnothing(plots) && push!(plots, p)
    end

    return container
end

function render!(
    container,
    protein::AbstractVector{Protein.Chain};
    colorscheme::ColorScheme = default_colorscheme,
    color_vectors::AbstractVector{<:AbstractVector{<:RGB}} = [colorscheme[LinRange(0, 1, length(chain))] for chain in protein],
    missing_residue_color = colorant"gray",
    kwargs...
)
    @assert Protein.has_assigned_ss(protein) "Protein must have assigned secondary structure."
    @assert length(protein) == length(color_vectors)
    @assert length.(protein) == length.(color_vectors)
    for (chain, colors) in zip(protein, color_vectors)
        subchain_ranges = split_by_resnum(chain)
        if length(subchain_ranges) == 1
            render!(container, chain, colors; kwargs...)
        else
            subchains = [chain[r] for r in subchain_ranges]
            for (subchain, r) in zip(subchains, subchain_ranges)
                render!(container, subchain, colors[r]; kwargs...)
            end
            draw_lines_between_subchains!(container, subchains, missing_residue_color; kwargs...)
        end
    end

    return container
end