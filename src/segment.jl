function segment_startpoint(segment::Segment)
    return segment[:,1,1] # first nitrogen atom
end

function segment_endpoint(segment::Segment)
    segment_stop = segment.range.stop
    if segment_stop == length(segment.chain)
        return segment.chain[:,3,segment_stop] # last carbon atom
    else
        return segment.chain[:,1,segment_stop+1] # next nitrogen atom
    end
end