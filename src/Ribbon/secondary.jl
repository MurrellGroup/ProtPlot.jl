import AssigningSecondaryStructure as ASS

const COIL = 1
const HELIX = 2
const STRAND = 3

const MIN_HELIX_LENGTH = 4
const MIN_STRAND_LENGTH = 2

# TODO: fix rendering edge cases so this can be removed
function clean_secondary_structure!(secondary_structure::Vector{Int})
    n = length(secondary_structure)
    i = 1
    while i <= n
        current_structure = secondary_structure[i]
        segment_start = i
        while i <= n && secondary_structure[i] == current_structure
            i += 1
        end
        segment_end = i - 1
        segment_length = segment_end - segment_start + 1
        for (ss, min_len) in [(HELIX, MIN_HELIX_LENGTH), (STRAND, MIN_STRAND_LENGTH)]
            if current_structure == ss && segment_length < min_len
                for j in segment_start:segment_end
                    secondary_structure[j] = COIL
                end
            end
        end
    end
    return secondary_structure
end

function _assign_secondary_structure(chain_backbones::Vector{Array{T,3}}) where T<:Real
    secondary_structure_by_chain = ASS.assign_secondary_structure(chain_backbones)
    clean_secondary_structure!.(secondary_structure_by_chain)
    return secondary_structure_by_chain
end

function segments(secondary_structure::Vector{Int})
    segment_ranges = Tuple{Int,UnitRange{Int}}[]
    start_i = 1
    for (i, ss) in enumerate(secondary_structure)
        if ss != secondary_structure[start_i]
            if ss != COIL && secondary_structure[start_i] != COIL
                push!(segment_ranges, (COIL, i:i-1))
            end
            push!(segment_ranges, (secondary_structure[start_i], start_i:i-1))
            start_i = i
        end
    end
    push!(segment_ranges, (secondary_structure[start_i], start_i:length(secondary_structure)))
    return segment_ranges
end
