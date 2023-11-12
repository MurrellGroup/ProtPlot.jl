module ProtPlot

using Backboner
export Backboner

using GLMakie
using Colors, ColorSchemes
using LinearAlgebra

include("utils.jl")
include("segment.jl")
include("splines.jl")
include("shapes.jl")
include("render.jl")

end
