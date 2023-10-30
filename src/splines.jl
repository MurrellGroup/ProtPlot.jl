import Dierckx

function spline(points_matrix::AbstractMatrix{<:Real}; m::Integer=10, k::Integer=3, N=nothing)
    L = size(points_matrix, 2)
    N = isnothing(N) ? L*m : N
    spl = Dierckx.ParametricSpline(1:L, points_matrix, k=k)
    linrange = Dierckx.LinRange(1, L, N)
    points_matrix_fine = Dierckx.evaluate(spl, linrange)
    return points_matrix_fine
end

function spline(points::AbstractVector{<:AbstractVector{<:Real}}; kwargs...)
    points_matrix = hcat(points...)
    return spline(points_matrix; kwargs...)
end
