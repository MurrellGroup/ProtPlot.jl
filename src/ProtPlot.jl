module ProtPlot

using Backboner
using Makie

import Backboner.Protein: readpdb

export Backboner, readpdb

include("ribbon/ribbon.jl")
include("attention/attention.jl")
include("ramachandran.jl")

end
 