using AssigningSecondaryStructure

@enum SSClass Unassigned Coil Helix Strand

const COIL_CODES = Set(['C', 'T', 'S', '-', ' '])
const HELIX_CODES = Set(['G', 'H', 'I'])
const STRAND_CODES = Set(['E', 'B'])

const SS_CLASS_DICT = Dict{Char, SSClass}()

for (codes, structure) in zip([COIL_CODES, HELIX_CODES, STRAND_CODES], [Coil, Helix, Strand])
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

"""
    assign_secondary_structure!(protein)

Uses a simplified version of DSSP to fill the secondary structure vector of each chain with '-' (coil/loop), 'H' (helix), and 'E' (strand).
"""
function AssigningSecondaryStructure.assign_secondary_structure!(protein::Vector{Protein.Chain})
    ss_vectors = assign_secondary_structure(Protein.ncaco_coords.(protein))
    for (chain, ssvector) in zip(protein, ss_vectors)
        @assert length(chain.ssvector) == length(ssvector)
        chain.ssvector .= get.(Ref(INT_TO_SS_CODE), ssvector, '-')
        clean_secondary_structure!(chain.ssvector)
    end
    return protein
end

"""
    assign_secondary_structure(protein)

Returns a new protein with secondary structure assigned.
"""
function AssigningSecondaryStructure.assign_secondary_structure(protein::Vector{Protein.Chain})
    new_protein = deepcopy(protein)
    assign_secondary_structure!(new_protein)
    return new_protein
end