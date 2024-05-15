# This code doesn't really belong in this package, but was left here for convenience.

module Attention

using ..ProtPlot
using Backboner
using Makie

export PointAttention

"""
Holds L 3D points and an HxLxL attention intensity tensor.
The attention to the first i residues at time i is intensities[:, i, 1:i].
Each column of this matrix represents the attention of some head to the corresponding points[1:i].
"""
struct PointAttention{T}
    points::Matrix{T}
    intensities::Array{T, 3}

    function PointAttention(points::AbstractMatrix{<:Real}, intensities::AbstractArray{<:Real, 3})
        size(points, 1) == 3 || throw(ArgumentError("The points matrix must have 3 rows."))
        l, h = size(points, 2), size(intensities, 1)
        size(intensities) == (h, l, l) || throw(ArgumentError("The size of the intensities tensor must match the number of points."))
        T = promote_type(eltype(points), eltype(intensities))
        new{T}(points, intensities)
    end
end

function draw_lines_from_point!(container,
    point::P, other_points::AbstractVector{P}, linewidths::AbstractVector{<:Real};
    color = RGB(1, 1, 1), linewidth_factor = 3.0, plot_list = nothing, kwargs...
) where P
    length(other_points) == length(linewidths) || throw(ArgumentError("The number of linewidths must match the number of other points.")) 
    xs, ys, zs = [reduce(vcat, ([point[i], other_point[i]] for other_point in other_points); init=Vector{eltype(point)}()) for i in 1:3]
    p = linesegments!(container, xs, ys, zs; linewidth=linewidths .* linewidth_factor, color, transparency=true, kwargs...)
    !isnothing(plot_list) && push!(plot_list, p)
end

"""
Take H points, and an LxH attention intensity matrix.
For each column slice of the attention matrix, draw lines to the corresponding points with the intensities in the vector that exceed a certain threshold. 
The thickness of the lines should be proportional to the value in the attention matrix.
"""
function draw_attention!(container,
    point::P, other_points::AbstractVector{P}, intensity_matrix::AbstractMatrix{<:Real};
    colors, colormap=:jet, threshold::Real=1.0, colorrange=(0,1), kwargs...
) where P
    h, l = size(intensity_matrix)
    #println(h, " ", l, " ", length(other_points))
    l == length(other_points) || throw(ArgumentError("The number of points must match the number of rows in the attention matrix."))
    for i in 1:h
        intensity_vector = @view intensity_matrix[i, :]
        js = intensity_vector .> threshold # TODO: alternatively get the K largest values
        draw_lines_from_point!(container, point, @view(other_points[js]), @view(intensity_vector[js]); color=colors[i], colormap, colorrange, kwargs...)
    end
    return nothing
end

function draw_attention_slice!(container, i::Int, attention::PointAttention; kwargs...)
    draw_attention!(container, eachcol(attention.points)[i], eachcol(attention.points)[1:i], @view(attention.intensities[:, i, 1:i]); kwargs...)
    return nothing
end

include("animate.jl")

end