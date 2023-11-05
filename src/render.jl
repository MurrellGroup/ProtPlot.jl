export render, render!

function render!(scene::Scene, segment::Segment{Loop})
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    controls = @view alphacarbon_coord_matrix(segment)[:, 2:end]
    coords = hcat(startpoint, controls, endpoint)
    surface_vertices, color_matrix = tube(coords, 0.3, spline_quality=10, tube_quality=10, color_start=segment.range.start/length(segment.chain), color_end=segment.range.stop/length(segment.chain))
    surface!(scene, eachslice(surface_vertices, dims=1)..., color=color_matrix)
end

function render!(scene::Scene, segment::Segment{Helix})
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    controls = @view alphacarbon_coord_matrix(segment)[:, 2:end] # startpoint is first point instead. including first N *and* CA could mess with normals
    coords = hcat(startpoint, controls, endpoint)
    surface_vertices, color_matrix = tube(coords, 1, spline_quality=10, tube_quality=10, x_elongation=0.25, color_start=segment.range.start/length(segment.chain), color_end=segment.range.stop/length(segment.chain))
    surface!(scene, eachslice(surface_vertices, dims=1)..., color=color_matrix)
end

function render!(scene::Scene, segment::Segment{Strand})
    startpoint = segment_startpoint(segment)
    endpoint = segment_endpoint(segment)
    oxygen_coords_side1 = @view oxygen_coord_matrix(segment)[:, 1:2:end-1]
    oxygen_coords_side2 = @view oxygen_coord_matrix(segment)[:, 2:2:end]
    coords1 = hcat(startpoint, oxygen_coords_side1, endpoint)
    coords2 = hcat(startpoint, oxygen_coords_side2, endpoint)
    surface_vertices, color_matrix = sheet(coords1, coords2, spline_quality=5, thickness=0.5, width_pad=0.5, color_start=segment.range.start/length(segment.chain), color_end=segment.range.stop/length(segment.chain))
    surface!(scene, eachslice(surface_vertices, dims=1)..., color=color_matrix)
end

function render!(scene::Scene, chain::Chain{4})
    @assert !has_missing_ss(chain)
    for segment in segments(chain)
        render!(scene, segment)
    end
end

function render(backbone::Backbone{4})
    has_missing_ss(backbone) && assign_secondary_structure!(backbone)
    scene = Scene(backgroundcolor=:black)#, axis=(;type=Axis3, aspect=:data))
    cam3d!(scene)
    for chain in backbone
        render!(scene, chain)
    end
    center!(scene)
    display(scene)
end

function backbone_gif(backbones::Vector{<:Backbone{4}})
    scene = Scene(backgroundcolor=:black)#, axis=(;type=Axis3, aspect=:data))
    record(scene, "backbone.gif", eachindex(backbones)) do i
        empty!(scene)
        cam3d!(scene)
        backbone = backbones[i]
        backbone[1].coords[1,:,:] += backbone[1].coords[1,:,:] * 0.001i 
        backbone[2].coords[1,:,:] += backbone[1].coords[1,:,:] * 0.001i 
        has_missing_ss(backbone) && assign_secondary_structure!(backbone)
        for chain in backbone
            render!(scene, chain)
        end
        update_cam!(scene)
    end
end