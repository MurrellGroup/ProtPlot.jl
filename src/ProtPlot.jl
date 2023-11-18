module ProtPlot

import Backboner: Backboner, pdb_to_protein
export Backboner, pdb_to_protein

include("ribbon/ribbon.jl")

import .Ribbon: ribbon, ribbon!
export Ribbon, ribbon, ribbon!

end
