include("shapes/shapes.jl")

function render!(
    ribbon::Ribbon,
    ::Val{:coil},
    segment_range::UnitRange{Int},
    chain::Protein.Chain,
    colors;
    radius = 0.2,
    spline_quality = 20,
    slice_quality = 20,
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
    surface!(ribbon, eachslice(surface_vertices, dims=1)..., color=color_matrix)
    nothing
end

function render!(
    ribbon::Ribbon,
    ::Val{:helix},
    segment_range::UnitRange{Int},
    chain::Protein.Chain,
    colors;
    helix_radius = 1.0,
    helix_width = 1.0,
    helix_thickness = 0.25,
    spline_quality = 20,
    slice_quality = 20,
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
    surface!(ribbon, eachslice(surface_vertices, dims=1)..., color=color_matrix)
    nothing
end

function render!(
    ribbon::Ribbon,
    ::Val{:strand},
    segment_range::UnitRange{Int},
    chain::Protein.Chain,
    colors;
    strand_width = 2.0,
    strand_thickness = 0.5,
    spline_quality = 20,
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
    surface!(ribbon, eachslice(surface_vertices, dims=1)..., color=color_matrix)
    nothing
end

function render!(
    ribbon::Ribbon,
    chain::Protein.Chain,
    colors::AbstractVector{<:RGB};
    kwargs...
)
    for (ss_name, segment_range) in segments(chain)
        render!(ribbon, Val(ss_name), segment_range, chain, colors; kwargs...)
    end
    nothing
end

function draw_lines_between_subchains!(ribbon::Ribbon, subchains::AbstractVector{Protein.Chain}, color::RGB; linewidth=2, kwargs...)
    for (i, j) in zip(1:length(subchains)-1, 2:length(subchains))
        startpoint, endpoint = subchains[i].backbone[end], subchains[j].backbone[begin]
        n_segments = trunc(Int, norm(endpoint - startpoint) / 0.8)
        xs, ys, zs = [LinRange(startpoint[i], endpoint[i], 2*n_segments) for i in 1:3]
        linesegments!(ribbon, xs, ys, zs; linewidth=linewidth, color=color, transparency=true)
    end
    nothing
end

function render!(ribbon::Ribbon, chains::Vector{Protein.Chain})
    for (chain, colors) in zip(chains, ribbon.attributes.color_vectors)
        subchain_ranges = split_by_resnum(chain)
        if length(subchain_ranges) == 1
            render!(ribbon, chain, colors; kwargs...)
        else
            subchains = [chain[r] for r in subchain_ranges]
            for (subchain, r) in zip(subchains, subchain_ranges)
                render!(ribbon, subchain, colors[r]; kwargs...)
            end
            draw_lines_between_subchains!(ribbon, subchains, missing_residue_color; kwargs...)
        end
    end

    return nothing
end