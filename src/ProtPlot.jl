module ProtPlot

using BackBoner
export BackBoner

import AssigningSecondaryStructure as ASS
export ASS

using GLMakie
using Colors, ColorSchemes
using LinearAlgebra

include("segment.jl")
include("splines.jl")
include("shapes.jl")
include("render.jl")

end
