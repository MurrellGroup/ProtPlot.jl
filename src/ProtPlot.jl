module ProtPlot

export Ribbon, ribbon, ribbon!
export ribbon_scene

export Ramachandran, ramachandran, ramachandran!

export readchains, pdbentry, @pdb_str

using Makie
using ColorTypes

import ProteinChains: ProteinChain, readchains, pdbentry, @pdb_str
import Backboner: Backbone, get_torsion_angles

include("Ribbon/Ribbon.jl")
include("ramachandran.jl")

end
