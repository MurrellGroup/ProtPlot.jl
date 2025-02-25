module ProtPlot

using Makie
using ColorTypes

using ProteinChains

import BioStructures

export pdbentry, @pdb_str

include("ribbon/ribbon.jl")
export Ribbon, ribbon, ribbon!
export ribbon_scene

include("atomplot.jl")
export AtomPlot, atomplot, atomplot!

include("ramachandran.jl")
export Ramachandran, ramachandran, ramachandran!

const _ProteinPlot = Union{Ribbon, AtomPlot, Ramachandran}

Makie.convert_arguments(P::Type{<:_ProteinPlot}, path::AbstractString) = Makie.convert_arguments(P, BioStructures.retrievepdb(path))

end
