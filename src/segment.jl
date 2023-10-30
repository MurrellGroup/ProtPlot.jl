struct Segment{C}
    range::UnitRange{Int}
    backbone::Backbone
    subbackbone::Backbone

    function Segment{C}(range::UnitRange{Int}, backbone::Backbone) where C
        return new{C}(range, backbone, backbone[range])
    end
end

function segmenter(ss::AbstractVector{ASS.SSClass}, backbone::Backbone)
    start_idx = 1
    end_idx = 1
    segments = Segment[]
    for (i, class) in enumerate(ss)
        if class != ss[start_idx]
            push!(segments, Segment{ss[start_idx]}(start_idx:end_idx, backbone))
            start_idx = i
        end
        end_idx = i
    end
    push!(segments, Segment{ss[start_idx]}(start_idx:end_idx, backbone))
    return segments
end

function segment_startpoint(segment::Segment)
    return segment.subbackbone[1][:,1] # first nitrogen atom
end

function segment_stoppoint(segment::Segment)
    stop = segment.range.stop
    if stop == length(segment.backbone)
        return segment.backbone[stop][:,3] # last carbon atom
    else
        return segment.backbone[stop+1][:,1] # next nitrogen atom
    end
end