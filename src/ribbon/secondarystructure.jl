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

function clean_secondary_structure!(ss_codes::Vector{Char})
    n = length(ss_codes)
    i = 1

    while i <= n
        current_structure = ss_codes[i]
        segment_start = i

        while i <= n && ss_codes[i] == current_structure
            i += 1
        end
        segment_end = i - 1
        segment_length = segment_end - segment_start + 1

        for (code, min_len) in [('H', MIN_HELIX_LENGTH), ('E', MIN_STRAND_LENGTH)]
            if current_structure == code && segment_length < min_len
                for j in segment_start:segment_end
                    ss_codes[j] = '-'
                end
            end
        end
    end

    return ss_codes
end

const INT_TO_SS_CODE = ['-', 'H', 'E']

function _assign_secondary_structure(chains::Vector{Protein.Chain})
    ss_numbers_by_chain = ASS.assign_secondary_structure(chains)
    ss_codes_by_chain = Vector{Char}[]
    for ss_numbers in ss_numbers_by_chain
        ss_codes = get.(Ref(INT_TO_SS_CODE), ss_numbers, '-')
        clean_secondary_structure!(ss_codes)
        push!(ss_codes_by_chain, ss_codes)
    end
    return ss_codes_by_chain
end

function _assign_secondary_structure!(chains::Vector{Protein.Chain})
    ss_codes_by_chain = _assign_secondary_structure(chains)
    for (chain, ss_codes) in zip(chains, ss_codes)
        chain.ssvector .= ss_codes
    end
    return chains
end

function segments(ss_codes::Vector{Char})
    ss_names = [SS_NAME_DICT[ss] for ss in ss_codes]
    segment_ranges = Tuple{Symbol, UnitRange{Int}}[]
    start_i = 1

    for (i, ss_name) in enumerate(ss_names)
        if ss_name != ss_names[start_i]
            if ss_name != :Coil && ss_names[start_i] != :Coil
                push!(segment_ranges, (:Coil, i:i-1))
            end
            push!(segment_ranges, (ss_names[start_i], start_i:i-1))
            start_i = i
        end
    end
    push!(segment_ranges, (ss_names[start_i], start_i:length(ss_names)))

    return segment_ranges
end

segments(chain::Protein.Chain) = segments(chain.ssvector)