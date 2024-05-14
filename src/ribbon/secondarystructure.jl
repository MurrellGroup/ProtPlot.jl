import AssigningSecondaryStructure as ASS

const COIL_CODES = Set(['C', 'T', 'S', '-', ' '])
const HELIX_CODES = Set(['G', 'H', 'I'])
const STRAND_CODES = Set(['E', 'B'])

const SS_NAME_DICT = Dict{Char, Symbol}()

for (codes, structure) in zip([COIL_CODES, HELIX_CODES, STRAND_CODES], [:coil, :helix, :strand])
    for code in codes
        SS_CLASS_DICT[code] = structure
    end
end

const MIN_HELIX_LENGTH = 4
const MIN_STRAND_LENGTH = 2

function clean_secondary_structure!(ss_vector::Vector{Char})
    n = length(ss_vector)
    i = 1

    while i <= n
        current_structure = ss_vector[i]
        start = i

        while i <= n && ss_vector[i] == current_structure
            i += 1
        end
        segment_end = i - 1
        segment_length = segment_end - start + 1

        for (code, max_len) in [('H', MIN_HELIX_LENGTH), ('E', MIN_STRAND_LENGTH)]
            if current_structure == code && segment_length < max_len
                for j in start:segment_end
                    ss_vector[j] = '-'
                end
            end
        end
    end

    return ss_vector
end

const INT_TO_SS_CODE = ['-', 'H', 'E']

function _assign_secondary_structure(chains::Vector{Protein.Chain})
    ssvectors = ASS.assign_secondary_structure(chains)
    for ssvector in ssvectors
        ssvector .= get.(Ref(INT_TO_SS_CODE), ssvector, '-')
        clean_secondary_structure!(ssvector)
    end
    return ssvectors
end

function _assign_secondary_structure!(chains::Vector{Protein.Chain})
    ssvectors = _assign_secondary_structure(chains)
    for (chain, ssvector) in zip(chains, ssvectors)
        chain.ssvector .= ssvector
    end
    return chains
end