import Dierckx

# TODO: look at residues next to the segment to get a smoother spline at the ends

function spline(points_matrix::AbstractMatrix{<:Real}; m::Integer=10, k::Integer=3, N=nothing)
    L = size(points_matrix, 2)
    N = isnothing(N) ? L*m : N
    spl = Dierckx.ParametricSpline(1:L, points_matrix, k=k)
    linrange = Dierckx.LinRange(1, L, N)
    points_matrix_fine = Dierckx.evaluate(spl, linrange)
    return points_matrix_fine
end

function spline(
    points_matrix::AbstractMatrix{<:Real}, ghost_control_start::AbstractVector{<:Real}, ghost_control_end::AbstractVector{<:Real};
    m::Integer=10, k::Integer=3, N=nothing
)
    L = size(points_matrix, 2)
    N = isnothing(N) ? L*m : N
    ghost_points_matrix = [ghost_control_start points_matrix ghost_control_end]
    spl = Dierckx.ParametricSpline(0:L+1, ghost_points_matrix, k=k)
    linrange = Dierckx.LinRange(1, L, N)
    points_matrix_fine = Dierckx.evaluate(spl, linrange)
    return points_matrix_fine
end

function spline(points::AbstractVector{<:AbstractVector{<:Real}}; kwargs...)
    points_matrix = hcat(points...)
    return spline(points_matrix; kwargs...)
end
