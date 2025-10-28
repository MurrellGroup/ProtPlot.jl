using ProteinChains: get_atoms, atom_coords, atom_symbol
using NearestNeighbors
using StaticArrays

const ATOM_COLORS = Dict(
    "C" => :gray60,
    "H" => :white,
    "O" => :red,
    "N" => :royalblue2,
    "S" => :yellow,
    "P" => :green,
    "X" => :darkred,
)

const ATOM_SIZE_FACTORS = Dict(
    "H" => 0.75f0,
)

const BOND_THRESHOLDS = Dict(
    # C–C bonds
    ("C", "C") => 1.54,
 
    # C–H bonds
    ("C", "H") => 1.09,
 
    # C–N bonds
    ("C", "N") => 1.47,
 
    # C–O bonds
    ("C", "O") => 1.43,
 
    # C–S bonds
    ("C", "S") => 1.82,
 
    # C–halogen bonds
    ("C", "F") => 1.35,
    ("C", "Cl") => 1.77,
    ("C", "Br") => 1.94,
    ("C", "I") => 2.14,
 
    # N–H and O–H
    ("N", "H") => 1.01,
    ("O", "H") => 0.96,
 
    # N–O and N–N
    ("N", "O") => 1.40,
    ("N", "N") => 1.45,
 
    # O–O and S–S
    ("O", "O") => 1.48,
    ("S", "S") => 2.05,
    ("S", "O") => 1.63,
 
    # Metal–ligand bonds (approximate, vary by oxidation state and coordination)
    ("Fe", "N") => 2.00,
    ("Fe", "O") => 1.65,
    ("Fe", "S") => 2.25,
    ("Cu", "N") => 2.00,
    ("Cu", "O") => 1.95,
    ("Cu", "S") => 2.25,
    ("Ni", "N") => 1.90,
    ("Ni", "O") => 1.85,
    ("Ni", "S") => 2.13,
    ("Co", "N") => 1.96,
    ("Co", "O") => 1.88,
    ("Zn", "N") => 2.05,
    ("Zn", "O") => 1.95,
    ("Mn", "O") => 1.90,
    ("Cr", "O") => 1.60,
    ("V", "O") => 1.58,
    ("Mo", "O") => 1.70,
    ("W", "O") => 1.71,
    ("Pt", "Cl") => 2.30,
    ("Pd", "Cl") => 2.28,
    ("Ag", "N") => 2.20,
    ("Ag", "O") => 2.10
)

get_bond_threshold(a, b, fallback=0.0) = max(get(BOND_THRESHOLDS, (a, b), fallback), get(BOND_THRESHOLDS, (b, a), fallback))

@recipe(AtomPlot, atoms) do scene
    Attributes(
        colormap = :jet,
        colorrange = nothing,
        color = ATOM_COLORS,
        default_color = :gray40,
        size_factor = ATOM_SIZE_FACTORS,
        default_size = 0.7f0,
        bond_width = 0.25f0,
        show_bonds = false,
        bond_threshold_tolerance = 0.0f0,
        alpha = 1f0,
        transform_marker = true,
    )
end

function get_from_dict(dict, atoms, default)
    dict isa Dict ? [get(dict, atom_symbol(atom), default) for atom in atoms] : dict
end

function get_pairs(positions, atoms, tolerance)
    tree = KDTree(positions)
    indices, distances = knn(tree, positions, min(5, length(positions)))
    pairs = Tuple{Int,Int}[]
    for i in eachindex(positions)
        for (j, d) in zip(indices[i], distances[i])
            if i < j && d < get_bond_threshold(atom_symbol(atoms[i]), atom_symbol(atoms[j])) + tolerance
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
    markersize = [Vec3f(1.0,1.0,x) for x in lengths]
    rotation = [Makie.rotation_between(SVector(0.0,0.0,1.0), d) for d in dirs]
    i_midpoints = (midpoints + p1) / 2
    j_midpoints = (midpoints + p2) / 2
    return (; midpoints, markersize, rotation, i_midpoints, j_midpoints)
end

function Makie.plot!(plot::AtomPlot{<:Tuple{<:AbstractVector{<:Atom}}})
    positions = @lift atom_coords.($(plot.atoms))
    markersize = @lift get_from_dict($(plot.size_factor), $(plot.atoms), 1.0f0) * $(plot.default_size)
    color = @lift get_from_dict($(plot.color), $(plot.atoms), $(plot.default_color))

    meshscatter!(plot, positions;
        markersize, color,
        plot.colormap,
        plot.colorrange,
        plot.alpha,
        plot.transform_marker,
    )

    pairs = @lift get_pairs(atom_coords.($(plot.atoms)), $(plot.atoms), $(plot.bond_threshold_tolerance))
    bonds = @lift get_bonds(atom_coords.($(plot.atoms)), $pairs, $(plot.show_bonds))

    !plot.show_bonds[] && return plot

    # need to recompute because {i,j}_color depending on pairs and color will fail
    # because the inputs need to be updated simultaneously. see 
    bond_colors = @lift begin
        atoms = $(plot.atoms)
        color = get_from_dict($(plot.color), atoms, $(plot.default_color))
        positions = atom_coords.(atoms)
        pairs = get_pairs(positions, atoms, $(plot.bond_threshold_tolerance))
        i_color = color isa AbstractVector ? [color[i] for (i, _) in pairs] : color
        j_color = color isa AbstractVector ? [color[j] for (_, j) in pairs] : color
        return (; i_color, j_color)
    end
    i_color = @lift $bond_colors.i_color
    j_color = @lift $bond_colors.j_color

    meshscatter!(plot, (@lift $bonds.i_midpoints);
        rotation = (@lift $bonds.rotation),
        marker = (@lift Cylinder(Point3f(0.0,0.0,-0.25), Point3f(0.0,0.0,0.25), $(plot.bond_width))),
        markersize = (@lift $bonds.markersize),
        color = i_color,
        #colormap = plot.colormap,
        plot.alpha, plot.transform_marker
    )

    meshscatter!(plot, (@lift $bonds.j_midpoints);
        rotation = (@lift $bonds.rotation),
        marker = (@lift Cylinder(Point3f(0.0,0.0,-0.25), Point3f(0.0,0.0,0.25), $(plot.bond_width))),
        markersize = (@lift $bonds.markersize),
        #color = j_color,
        #colormap = plot.colormap,
        plot.alpha, plot.transform_marker
    )

    return plot
end

Makie.args_preferred_axis(::Type{<:AtomPlot}, atoms) = LScene

Makie.convert_arguments(::Type{<:AtomPlot}, atoms::AbstractVector{<:AbstractVector{<:Atom}}) = (reduce(vcat, atoms),)
Makie.convert_arguments(::Type{<:AtomPlot}, chain::Union{ProteinChain,Backbone}) = Makie.convert_arguments(AtomPlot, get_atoms(chain))
Makie.convert_arguments(::Type{<:AtomPlot}, struc::ProteinStructure) = Makie.convert_arguments(AtomPlot, vcat(struc.atoms, reduce(vcat, get_atoms.(struc))...))
Makie.convert_arguments(::Type{<:AtomPlot}, frames::Frames) = Makie.convert_arguments(AtomPlot, Backbone(frames(STANDARD_RESIDUE)))

Makie.convert_arguments(::Type{<:AtomPlot}, struc::BioStructures.StructuralElementOrList) = (convert.(Atom{Float64}, BioStructures.collectatoms(struc)),)