function render!(segment::Segment{ASS.Loop})
    startpoint = segment_startpoint(segment)
    stoppoint = segment_stoppoint(segment)
    anchors = alphacarbon_coord_matrix(segment.subbackbone[2:end])
    coords = hcat(startpoint, anchors, stoppoint)
    surface_vertices, color_matrix = tube(coords, 0.3, spline_quality=10, tube_quality=10, color_start=segment.range.start/length(segment.backbone), color_end=segment.range.stop/length(segment.backbone))
    surface!(eachslice(surface_vertices, dims=1)..., color=color_matrix)
end

function render!(segment::Segment{ASS.Helix})
    startpoint = segment_startpoint(segment)
    stoppoint = segment_stoppoint(segment)
    anchors = alphacarbon_coord_matrix(segment.subbackbone[2:end]) # startpoint is first point instead. including first N *and* CA could mess with normals
    coords = hcat(startpoint, anchors, stoppoint)
    surface_vertices, color_matrix = tube(coords, 1, x_elongation=0.15, color_start=segment.range.start/length(segment.backbone), color_end=segment.range.stop/length(segment.backbone))
    surface!(eachslice(surface_vertices, dims=1)..., color=color_matrix)
end

function render!(segment::Segment{ASS.Strand})
    startpoint = segment_startpoint(segment)
    stoppoint = segment_stoppoint(segment)
    oxygen_coords = oxygen_coord_matrix(segment.subbackbone)
    oxygen_coords_side1 = oxygen_coords[:, 1:2:end-1]
    oxygen_coords_side2 = oxygen_coords[:, 2:2:end]
    coords1 = hcat(startpoint, oxygen_coords_side1, stoppoint)
    coords2 = hcat(startpoint, oxygen_coords_side2, stoppoint)
    surface_vertices, color_matrix = sheet(coords1, coords2, thickness=0.5, width_pad=0.5, color_start=segment.range.start/length(segment.backbone), color_end=segment.range.stop/length(segment.backbone))
    surface!(eachslice(surface_vertices, dims=1)..., color=color_matrix)
end

function render!(backbone::Backbone{4}, ssclasses::Vector{ASS.SSClass})
    segments = segmenter(ssclasses, backbone)
    render!.(segments)
end

function render!(backbones::Vector{<:Backbone{4}})
    set_theme!(backgroundcolor = :black)

    ssclasses_vec = ASS.dssp([bb.coords for bb in backbones])
    scene = lines(Float64[])#, axis=(;type=Axis3, aspect=:data))
    for (backbone, ssclasses) in zip(backbones, ssclasses_vec)
        render!(backbone, ssclasses)
    end
    display(scene)
end