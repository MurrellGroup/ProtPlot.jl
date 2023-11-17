export Segment, extend_segment, segments

"""
    Segment{SS, T}

A segment of a chain with uniform secondary structure.
The segment can have mixed secondary structure if the `SS` type parameter is `Unassigned`.
"""
struct Segment{SS, T}
    chain::Chain{T}
    range::UnitRange{Int}
    backbone::Backbone{4,T}

    function Segment{SS}(chain::Chain{T}, range::UnitRange{Int}) where {SS, T}
        ssvector = view(chain.ssvector, range)
        @assert SS == Unassigned || all(==(SS), ssvector) "All residues in the '$SS' segment must have the '$SS' secondary structure"
        backbone = chain.backbone[range]
        return new{SS, T}(chain, range, backbone)
    end

    function Segment(chain::Chain{T}, range::UnitRange{Int}) where T
        SS = allequal(ssvector) ? ssvector[1] : Unassigned
        return Segment{SS}(chain, range)
    end
end

@inline Base.length(segment::Segment) = size(segment.backbone, 3)
@inline Base.size(segment::Segment) = (length(segment),)

Base.summary(segment::Segment{SS}) where SS = "$SS Segment of Chain $(segment.chain.id) with $(length(segment)) residues"
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
    #checkbounds(segment.chain.backbone, parent_vec_range)
    return Segment{Unassigned}(segment.chain, adjusted_range)
end

"""
    segments(chain)

Returns an array of segments of a chain.
The segments are defined by the secondary structure of the residues.
A chain with missing secondary structure information will throw an error.
"""
function segments(chain::Chain)
    has_missing_ss(chain) && @warn "Chain $(chain.id) has missing secondary structure information"

    ssvector = chain.ssvector
    start_idx = 1
    end_idx = 1
    segments = Segment[]
    for (i, ss) in enumerate(ssvector)
        if ss != ssvector[start_idx]
            push!(segments, Segment{ssvector[start_idx]}(chain, start_idx:end_idx))
            start_idx = i
        end
        end_idx = i
    end
    push!(segments, Segment{ssvector[start_idx]}(chain, start_idx:end_idx))
    return segments
end


# utilities

function segment_startpoint(segment::Segment, A=1) # nitrogen
    return segment.backbone.coords[:,A,1] # first residue
end

function segment_endpoint(segment::Segment, A=1) # nitrogen
    segment_stop = segment.range.stop
    if segment_stop == length(segment.chain)
        return segment.chain.backbone[:,3,segment_stop] # carbon of last residue
    else
        return segment.chain.backbone[:,A,segment_stop+1] # next residue
    end
end

function remove_singleton_strands!(chain::Chain)
    ssvector = chain.ssvector
    for i in 2:length(ssvector)-1
        if ssvector[i-1] != Strand && ssvector[i] == Strand && ssvector[i+1] != Strand
            ssvector[i] = Loop
        end
    end
end