@recipe(SpatialGraphPlot, positions, edges) do scene
    Attributes(
        linewidth = 1f0,
        alpha = 1f0,
        color = :violet,
    )
end

const Positions = AbstractVector{<:AbstractVector{<:Real}}
const Edges = AbstractVector{<:Tuple{Int,Int}}

function Makie.plot!(spatialgraph::SpatialGraphPlot{<:Tuple{Positions,Edges}})
    edge_positions = @lift map($(spatialgraph.edges)) do edge
        (Point3f($(spatialgraph.positions)[edge[1]]), Point3f($(spatialgraph.positions)[edge[2]]))
    end
    linesegments!(spatialgraph, spatialgraph.attributes, edge_positions;
        fxaa=true, transparency=true)
end

Makie.args_preferred_axis(::Type{<:SpatialGraphPlot}, positions, edges) = LScene

function Makie.convert_arguments(P::Type{<:SpatialGraphPlot}, chains::AbstractVector{<:ProteinChain}, g)
    positions = @views eachcol(hcat(map(c -> Backbone(c).coords[:,2:3:end], chains)...))
    Makie.convert_arguments(P, positions, g)
end

function Makie.convert_arguments(::Type{<:SpatialGraphPlot}, p::Positions, graph::AbstractMatrix{Bool})
    (p, Tuple.(findall(graph)))
end

import Graphs

function Makie.convert_arguments(::Type{<:SpatialGraphPlot}, p::Positions, graph::Graphs.AbstractGraph)
    (p, Tuple.(Graphs.edges(graph)))
end