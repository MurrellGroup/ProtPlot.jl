export Segment, extend_segment, segments

"""
    Segment{SS}

A segment of a chain with uniform secondary structure.
"""
struct Segment{SS}
    chain::Protein.Chain
    range::UnitRange{Int}
    backbone::Backbone

    function Segment{SS}(chain::Protein.Chain, range::UnitRange{Int}) where SS
        ssvector = view(chain.ssvector, range)
        ss_class_vector = [SS_CLASS_DICT[ss] for ss in ssvector]
        @assert SS == Unassigned || all(==(SS), ss_class_vector) "All residues in the '$SS' segment must have the '$SS' secondary structure"
        backbone = chain.backbone[3*range.start-2:3*range.stop]
        return new{SS}(chain, range, backbone)
    end

    function Segment(chain::Protein.Chain, range::UnitRange{Int})
        SS = SS_CLASS_DICT[chain.ssvector[range][1]]
        return Segment{SS}(chain, range)
    end
end

@inline Base.:(==)(segment1::Segment, segment2::Segment) = segment1.chain == segment2.chain && segment1.range == segment2.range
@inline Base.length(segment::Segment) = size(segment.backbone, 3)
@inline Base.size(segment::Segment) = (length(segment),)

Base.summary(segment::Segment{SS}) where SS = "$SS Segment of Protein.Chain $(segment.chain.id) with $(length(segment)) residues"
Base.show(io::IO, segment::Segment) = print(io, summary(segment))

@inline Base.getindex(segment::Segment{SS}, r::UnitRange{Int}) where SS = Segment{SS}(segment.chain, segment.range[r])

"""
    extend_segment(segment, range)

Returns a segment of the parent chain, extended to the given range.
If `segment` covers indices 3:4 of the parent chain, then `extend_segment(segment, 0:3)` will return a segment that covers indices 2:5, since 0 is one less than 1, and 3 is one more than 4.
`extend_segment(segment, 1:2)` would therefore return the same segment as `segment`.
This function is useful if one wishes to access the coordinates of the atoms of the parent chain that are adjacent to the segment.

!!! note
    The new segment will have missing secondary structure, since segments are meant to describe uniform secondary structure of contiguous residues.
"""
@inline function extend_segment(segment::Segment{SS}, range::UnitRange{Int}) where SS
    offset = segment.range.start - 1
    parent_vec_range = range .+ offset
    adjusted_range = max(1, parent_vec_range.start):min(length(segment.chain), parent_vec_range.stop)
    return Segment{Unassigned}(segment.chain, adjusted_range)
end

"""
    segments(chain)

Returns an array of segments of a chain.
The segments are defined by the secondary structure of the residues.
A chain with missing secondary structure information will throw an error.
"""
function segments(chain::Protein.Chain)
    !Protein.has_assigned_ss(chain) && @warn "Protein.Chain $(chain.id) has missing secondary structure information"

    ssvector = chain.ssvector
    ss_class_vector = [SS_CLASS_DICT[ss] for ss in ssvector]
    start_idx = 1
    end_idx = 1
    segments = Segment[]
    for (i, ss) in enumerate(ss_class_vector)
        if ss != ss_class_vector[start_idx]
            push!(segments, Segment{ss_class_vector[start_idx]}(chain, start_idx:end_idx))
            start_idx = i
        end
        end_idx = i
    end
    push!(segments, Segment{ss_class_vector[start_idx]}(chain, start_idx:end_idx))
    return segments
end


function segment_startpoint(segment::Segment)
    return segment.backbone[1] # nitrogen of first residue
end

function segment_endpoint(segment::Segment)
    if segment.range.stop == length(segment.chain)
        return segment.backbone[end] # carbonyl of last residue
    else
        return segment.chain.backbone[3*segment.range.stop+1] # next residue nitrogen
    end
end