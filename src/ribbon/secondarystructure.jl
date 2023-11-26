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