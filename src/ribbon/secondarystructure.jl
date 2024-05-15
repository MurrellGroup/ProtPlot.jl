import AssigningSecondaryStructure as ASS

const COIL_CODES = Set(['C', 'T', 'S', '-', ' '])
const HELIX_CODES = Set(['G', 'H', 'I'])
const STRAND_CODES = Set(['E', 'B'])

const SS_NAME_DICT = Dict{Char, Symbol}()

for (codes, structure) in zip([COIL_CODES, HELIX_CODES, STRAND_CODES], [:Coil, :Helix, :Strand])
    for code in codes
        SS_NAME_DICT[code] = structure
    end
end

const MIN_HELIX_LENGTH = 4
const MIN_STRAND_LENGTH = 2

function clean_secondary_structure!(ssvector::Vector{Char})
    n = length(ssvector)
    i = 1

    while i <= n
        current_structure = ssvector[i]
        segment_start = i

        while i <= n && ssvector[i] == current_structure
            i += 1
        end
        segment_end = i - 1
        segment_length = segment_end - segment_start + 1

        for (code, min_len) in [('H', MIN_HELIX_LENGTH), ('E', MIN_STRAND_LENGTH)]
            if current_structure == code && segment_length < min_len
                for j in segment_start:segment_end
                    ssvector[j] = '-'
                end
            end
        end
    end

    return ssvector
end

const INT_TO_SS_CODE = ['-', 'H', 'E']

function _assign_secondary_structure(chains::Vector{Protein.Chain})
    ssvectors_int = ASS.assign_secondary_structure(chains)
    ssvectors_char = Vector{Char}[]
    for ssvector_int in ssvectors_int
        ssvector_char = get.(Ref(INT_TO_SS_CODE), ssvector_int, '-')
        clean_secondary_structure!(ssvector_char)
        push!(ssvectors_char, ssvector_char)
    end
    return ssvectors_char
end

function _assign_secondary_structure!(chains::Vector{Protein.Chain})
    ssvectors = _assign_secondary_structure(chains)
    for (chain, ssvector) in zip(chains, ssvectors)
        chain.ssvector .= ssvector
    end
    return chains
end

function segments(chain::Protein.Chain)
    ssvector = chain.ssvector
    ss_names = [SS_NAME_DICT[ss] for ss in ssvector]
    segment_ranges = Tuple{Symbol, UnitRange{Int}}[]
    start_idx = 1

    for i in 2:length(ss_names)
        if ss_names[i] != ss_names[start_idx]
            push!(segment_ranges, (ss_names[start_idx], start_idx:i-1))
            start_idx = i
        end
    end

    # Push the last segment
    push!(segment_ranges, (ss_names[start_idx], start_idx:length(ss_names)))

    return segment_ranges
end