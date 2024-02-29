module ProtPlot

import Backboner
import Backboner.Protein: readpdb

export Backboner, readpdb

include("ribbon/ribbon.jl")

import .Ribbon: ribbon, ribbon!
export Ribbon, ribbon, ribbon!

end
