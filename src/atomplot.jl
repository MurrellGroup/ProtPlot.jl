using ProteinChains: get_atoms, atom_coords, atom_symbol
using NearestNeighbors
using StaticArrays

atom_colors = Dict(
    "C" => :gray60,
    "H" => :white,
    "O" => :red,
    "N" => :royalblue2,
    "S" => :yellow,
    "P" => :green,
)

@recipe(AtomPlot, atoms) do scene
    Attributes(
        colormap = :jet,
        colorrange = nothing,
        color = atom_colors,
        default_color = :gray40,
        size = Dict(),
        default_size = 0.7f0,
        bond_width = 0.25f0,
        bond_length_threshold = 2.0f0,
        show_bonds = false,
        alpha = 1f0,
        transform_marker = true,
    )
end

get_from_dict(dict, atoms, default) = dict isa Dict ? [get(dict, atom_symbol(atom), default) for atom in atoms] : dict

function get_pairs(positions, threshold)
    tree = KDTree(positions)
    indices, distances = knn(tree, positions, 5)
    pairs = Tuple{Int,Int}[]
    for i in eachindex(positions)
        for (j, d) in zip(indices[i], distances[i])
            if i < j && d < threshold
                push!(pairs, (i, j))
            end
        end
    end
    return pairs
end

function get_bonds(positions, pairs, show_bonds)
    pairs = show_bonds ? pairs : Tuple{Int,Int}[]
    p1 = @view positions[first.(pairs)]
    p2 = @view positions[last.(pairs)]
    midpoints = (p1 + p2) / 2
    displacements = p2 - p1
    dirs = normalize.(displacements)
    lengths = norm.(displacements)
    rotation = [Makie.rotation_between(SVector(0.0,0.0,1.0), d) for d in dirs]
    return (; midpoints, lengths, rotation)
end

function Makie.plot!(plot::AtomPlot{<:Tuple{<:AbstractVector{<:Atom}}})
    positions = @lift atom_coords.($(plot.atoms))
    markersize = @lift get_from_dict($(plot.size), $(plot.atoms), $(plot.default_size))
    color = @lift get_from_dict($(plot.color), $(plot.atoms), $(plot.default_color))

    meshscatter!(plot, positions;
        markersize, color,
        plot.colormap,
        plot.colorrange,
        plot.alpha,
        plot.transform_marker,
    )

    pairs = @lift get_pairs($positions, $(plot.bond_length_threshold))
    bonds = @lift get_bonds($positions, $pairs, $(plot.show_bonds))

    meshscatter!(plot, (@lift $bonds.midpoints);
        rotation = (@lift $bonds.rotation),
        markersize = (@lift $bonds.lengths),
        marker = (@lift Cylinder(Point3f(0.0,0.0,-0.5), Point3f(0.0,0.0,0.0), $(plot.bond_width))),
        color = (@lift $color isa AbstractVector ? [$color[i] for (i, _) in $pairs] : $color),
        plot.alpha,
        plot.transform_marker,
    )

    meshscatter!(plot, (@lift $bonds.midpoints);
        rotation = (@lift $bonds.rotation),
        markersize = (@lift $bonds.lengths),
        marker = (@lift Cylinder(Point3f(0.0,0.0,0.0), Point3f(0.0,0.0,0.5), $(plot.bond_width))),
        color = (@lift $color isa AbstractVector ? [$color[j] for (_, j) in $pairs] : $color),
        plot.alpha,
        plot.transform_marker,
    )

    return plot
end

Makie.args_preferred_axis(::Type{<:AtomPlot}, atoms) = LScene

Makie.convert_arguments(::Type{<:AtomPlot}, atoms::AbstractVector{<:AbstractVector{<:Atom}}) = (reduce(vcat, atoms),)
Makie.convert_arguments(::Type{<:AtomPlot}, chain::Union{ProteinChain,Backbone}) = Makie.convert_arguments(AtomPlot, get_atoms(chain))
Makie.convert_arguments(::Type{<:AtomPlot}, struc::ProteinStructure) = Makie.convert_arguments(AtomPlot, vcat(struc.atoms, reduce(vcat, get_atoms.(struc))...))
Makie.convert_arguments(::Type{<:AtomPlot}, frames::Frames) = Makie.convert_arguments(AtomPlot, Backbone(frames(STANDARD_RESIDUE)))

Makie.convert_arguments(::Type{<:AtomPlot}, struc::BioStructures.StructuralElementOrList) = (convert.(Atom{Float64}, BioStructures.collectatoms(struc)),)