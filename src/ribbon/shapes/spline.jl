import Dierckx

#="""
    spline(points_matrix)

Creates a spline curve from a matrix of points. The points are the columns of the matrix.
"""
function spline(points_matrix::AbstractMatrix{<:Real}; m::Integer=10, k::Integer=3, N=nothing)
    L = size(points_matrix, 2)
    N = isnothing(N) ? L*m : N
    spl = Dierckx.ParametricSpline(1:L, points_matrix, k=k)
    linrange = Dierckx.LinRange(1, L, N)
    points_matrix_fine = Dierckx.evaluate(spl, linrange)
    return points_matrix_fine
end=#

"""
    spline(points_matrix, ghost_control_start=nothing, ghost_control_end=nothing)

Allows for "ghost" control points at the start and end of the spline, to control the curvature at the ends.
"""
function spline(
    points_matrix::AbstractMatrix{<:Real},
    ghost_control_start::Union{Nothing, <:AbstractVector{<:Real}} = nothing,
    ghost_control_end::Union{Nothing, <:AbstractVector{<:Real}} = nothing;
    m::Integer=10, k::Integer=3, N=nothing
)
    has_start = !isnothing(ghost_control_start)
    has_end = !isnothing(ghost_control_end)
    L = size(points_matrix, 2)
    N = isnothing(N) ? L*m : N
    has_start && (points_matrix = hcat(ghost_control_start, points_matrix))
    has_end && (points_matrix = hcat(points_matrix, ghost_control_end))
    @assert size(points_matrix, 2) == length((has_start ? 0 : 1):L+(has_end ? 1 : 0))
    spl = Dierckx.ParametricSpline((has_start ? 0 : 1):L+(has_end ? 1 : 0), points_matrix, k=k)
    points_matrix_fine = Dierckx.evaluate(spl, LinRange(1, L, N))
    return points_matrix_fine
end

function spline(points::AbstractVector{<:AbstractVector{<:Real}}; kwargs...)
    points_matrix = hcat(points...)
    return spline(points_matrix; kwargs...)
end
