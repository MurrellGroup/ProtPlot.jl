module ProtPlot

export Ribbon, ribbon, ribbon!
export ribbon_scene

export Ramachandran, ramachandran, ramachandran!

export pdbentry, @pdb_str

using Makie
using ColorTypes

using ProteinChains
using Backboner: Backbone, get_torsion_angles

include("Ribbon/Ribbon.jl")
include("ramachandran.jl")

end
