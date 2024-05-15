module ProtPlot

using Backboner
using Makie

import Backboner.Protein: readpdb

export Backboner, readpdb

include("ribbon/ribbon.jl")
include("ramachandran.jl")

# TODO: move attention submodule elsewhere
include("attention/attention.jl")
using .Attention
export Attention, animate_attention

end
 