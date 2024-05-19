using LinearAlgebra

path_tangent(p1, p2, p3) = normalize!(p3 - p1)

import Dierckx

function spline(
    points::AbstractMatrix{<:Real};
    N::Integer = 10,
    k::Integer = min(3, size(points, 2)-1),
    r::UnitRange{Int} = 1:size(points, 2),
)
    L = size(points, 2)
    spl = Dierckx.ParametricSpline(1:L, points; k)
    return Dierckx.evaluate(spl, range(r.start, r.stop, N))
end

include("coil.jl")
include("helix.jl")
include("strand.jl")