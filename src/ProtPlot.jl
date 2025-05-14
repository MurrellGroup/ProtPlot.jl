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

include("spatialgraph.jl")
export SpatialGraphPlot, spatialgraphplot, spatialgraphplot!

include("trajectory.jl")
export animate_trajectory

const _ProteinPlot = Union{Ribbon, AtomPlot, Ramachandran, SpatialGraphPlot}

Makie.convert_arguments(P::Type{<:_ProteinPlot}, path::AbstractString, args...) = Makie.convert_arguments(P, BioStructures.retrievepdb(path), args...)

end
