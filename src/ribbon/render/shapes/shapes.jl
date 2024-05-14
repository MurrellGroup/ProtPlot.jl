path_tangent(p1, p2, p3) = normalize!(p3 - p1)
curved_path_normal(p1, p2, p3) = -normalize!(p1 - p2 + p3 - p2)

include("spline.jl")
include("coil.jl")
include("helix.jl")
include("strand.jl")