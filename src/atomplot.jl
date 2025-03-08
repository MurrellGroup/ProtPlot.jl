using ProteinChains: get_atoms, atom_coords, atom_symbol

atom_colors = Dict(
    "C" => :gray30,
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
        alpha=1f0,
    )
end

function Makie.plot!(atomplot::AtomPlot{<:Tuple{<:AbstractVector{<:Atom}}})
    positions = @lift atom_coords.($(atomplot.atoms))

    markersize = @lift if $(atomplot.size) isa Dict
        [get($(atomplot.size), atom_symbol(atom), $(atomplot.default_size)) for atom in $(atomplot.atoms)]
    else
        $(atomplot.size)
    end
 
    color = @lift if $(atomplot.color) isa Dict
        [get($(atomplot.color), atom_symbol(atom), $(atomplot.default_color)) for atom in $(atomplot.atoms)]
    else
        $(atomplot.color)
    end

    meshscatter!(atomplot, positions;
        markersize, color,
        atomplot.colormap,
        atomplot.colorrange,
        atomplot.alpha,
        transform_marker=true, # prevents constant apparant size in Axis3 containers
    )

    return atomplot
end

Makie.args_preferred_axis(::Type{<:AtomPlot}, atoms) = LScene

Makie.convert_arguments(::Type{<:AtomPlot}, atoms::AbstractVector{<:AbstractVector{<:Atom}}) = (reduce(vcat, atoms),)
Makie.convert_arguments(::Type{<:AtomPlot}, chain::Union{ProteinChain,Backbone}) = Makie.convert_arguments(AtomPlot, get_atoms(chain))
Makie.convert_arguments(::Type{<:AtomPlot}, struc::ProteinStructure) = Makie.convert_arguments(AtomPlot, vcat(struc.atoms, reduce(vcat, get_atoms.(struc))...))
Makie.convert_arguments(::Type{<:AtomPlot}, frames::Frames) = Makie.convert_arguments(AtomPlot, Backbone(frames(STANDARD_RESIDUE)))

Makie.convert_arguments(::Type{<:AtomPlot}, struc::BioStructures.StructuralElementOrList) = (convert.(Atom{Float64}, BioStructures.collectatoms(struc)),)