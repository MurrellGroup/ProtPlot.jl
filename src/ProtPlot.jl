module ProtPlot

using GLMakie
using ColorSchemes

import Backboner
import Backboner.Protein: readpdb

export Backboner, readpdb

include("ribbon/ribbon.jl")

using .Ribbon
export Ribbon, ribbon, ribbon!

include("attention/attention.jl")

using .Attention
export Attention, PointAttention, AttentionPlotIterator

include("animate.jl")

using .Animate
export Animate, animate_attention

include("ramachandran.jl")

end
 