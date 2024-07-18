module ProtPlot

using Makie
using ColorTypes

using Backboner
import Backboner.Protein: readpdb

export Backboner, readpdb

include("Ribbon/Ribbon.jl")
include("ramachandran.jl")

end
