function wider(heights)
    x = repeat(1:size(heights, 2), inner=size(heights, 1))
    group = repeat(1:size(heights, 1), outer=size(heights, 2))
    heights = reshape(heights, :)
    return x, group, heights
end
indexed_frames(locs, rots, step, inds) = Frames(rots[step][:,:,inds], 10 .* locs[step][:,inds])
triplicate(x) = repeat(x, inner = 3)
function animate_trajectory(export_path, samp::ProtPlot.ProteinStructure, trajectory;
                            aa_inds = trues(sum(length, samp)), pos_inds = trues(sum(length, samp)), 
                            rotation = 0.03, end_rotation_speedup = 0, size = (1280, 720), framerate = 22,
                            theme = :black, color_by_chain = false, atom_colormap = :jet,
                            kwargs...)
    ts, xt_locs,xt_rots,xt_aas,x̂1_locs,x̂1_rots,x̂1_aas = trajectory
    chainids = reduce(vcat, fill(i, length(chain)) for (i, chain) in enumerate(samp))
    steps = length(ts)
    frozen_prot = masked_out_structure(samp, pos_inds)
    theme == :black && set_theme!(theme_black())
    step = Observable(Int(round(steps/2)))
    timestep = Observable(1)
    fig = Makie.Figure(; size, figure_padding = 1)
    xtax = Axis3(fig[1:5, 1], perspectiveness=0.2, protrusions = (0,0,0,0), aspect=:data, width = 1000, height = 1000, tellwidth = false, tellheight = false; kwargs...)
    x̂1ax = Axis3(fig[1:5, 2], perspectiveness=0.2, protrusions = (0,0,0,0), aspect=:data, width = 1000, height = 1000, tellwidth = false, tellheight = false; kwargs...)

    Label(fig[1:5, 1, Top()], "xₜ", valign = :bottom, font = :bold, fontsize = 20, padding = (0, 0, 0, 0))
    Label(fig[1:5, 2, Top()], "x̂₁", valign = :bottom, font = :bold, fontsize = 20, padding = (0, 0, 0, 0))
    
    AAax = Axis(fig[6, 1:2])
    hidespines!(xtax)
    hidespines!(x̂1ax)
    hidedecorations!(xtax)
    hidedecorations!(x̂1ax)
    xtax.azimuth[] = -0.2 #pi/2
    xtax.elevation[] = 0.0 #pi/2
    
    x̂1ax.azimuth[] = -0.2 #pi/2
    x̂1ax.elevation[] = 0.0 #pi/2
    colorₜ = if color_by_chain
        triplicate((chainids .- 1) ./ (length(samp) - 1))
    else
        @lift triplicate(xt_aas[$step][aa_inds]) ./ 21
    end
    #Xt frames:
    xtframes0 = indexed_frames(xt_locs, xt_rots, 1, pos_inds)
    xtframesₜ = @lift indexed_frames(xt_locs, xt_rots, $step, pos_inds)
    fixed_xtribbon = ((length(frozen_prot) > 0) && (sum(length.(frozen_prot)) > 0)) ? ribbon!(xtax.scene, frozen_prot, colormap = :matter) : plot!(zeros(0))
    xtribbon = ribbon!(xtax, samp) #We need these to define the full expanse of the final state
    xtp = atomplot!(xtax.scene, xtframesₜ, color = colorₜ, colormap=atom_colormap, colorrange = (0.0, 1.0));
    #X̂1 frames:
    x̂1frames0 = indexed_frames(x̂1_locs, x̂1_rots, 1, pos_inds)
    x̂1framesₜ = @lift indexed_frames(x̂1_locs, x̂1_rots, $step, pos_inds)
    fixed_x̂1ribbon = ((length(frozen_prot) > 0) && (sum(length.(frozen_prot)) > 0)) ? ribbon!(x̂1ax.scene, frozen_prot, colormap = :matter) : plot!(zeros(0))
    x̂1ribbon = ribbon!(x̂1ax, samp) #We need these to define the full expanse of the final state
    x̂1p = atomplot!(x̂1ax.scene, x̂1framesₜ, color = colorₜ, colormap=atom_colormap, colorrange = (0.0, 1.0));
    #AA bar plot:
    aatitle = @lift $timestep <= length(ts) ? "t = $(rpad(round(ts[$timestep], digits = 3), 5, "0"))" : "t = $(rpad(round(1.0, digits = 3), 5, "0"))"
    Label(fig[6, 1:2, Top()], aatitle, valign = :bottom, font = :bold, padding = (0, 0, 0, 0))
    x,grp,freq = wider(x̂1_aas[1][:,aa_inds])
    freq_obs = @lift (x̂1_aas[$step][:,aa_inds])[:]
    Makie.barplot!(AAax, x, freq_obs, stack = grp, color = [2, 18, 1, 17, 12, 15, 13, 4, 16, 3, 5, 8, 10, 20, 14, 9, 6, 11, 19, 7, 21][grp], colormap = :turbo)
    colgap!(fig.layout, 0)
    rowgap!(fig.layout, 0)
    record(fig, export_path, -5:length(ts)+150, framerate=framerate) do t
        if t == -5
            step[] = 1
            timestep[] = 1
            delete!(xtax, xtribbon)
            delete!(x̂1ax, x̂1ribbon)
        end
        if 1 <= t <= length(ts)
            step[] = t
            timestep[] = t
        end
        if t == length(ts) + 1
            timestep[] = t
        end
        if t == length(ts) + 20
            delete!(xtax, fixed_xtribbon)
            delete!(x̂1ax, fixed_x̂1ribbon)
            ribbon!(xtax.scene, samp)
            ribbon!(x̂1ax.scene, samp, colors=[repeat([i/length(samp)], length(samp[i])) for i in 1:length(samp)], colormap=:hsv)
            ((length(frozen_prot) > 0) && (sum(length.(frozen_prot)) > 0)) ? ribbon!(xtax.scene, frozen_prot, colormap = :matter) : plot!(zeros(0))
            ((length(frozen_prot) > 0) && (sum(length.(frozen_prot)) > 0)) ? ribbon!(x̂1ax.scene, frozen_prot, colormap = :matter) : plot!(zeros(0))
        end
        if t == length(ts)+30
            delete!(xtax, xtp)
            delete!(x̂1ax, x̂1p)
        end
        x̂1ax.azimuth[] += rotation
        xtax.azimuth[] += rotation
        if t >= length(ts)+50
            x̂1ax.azimuth[] += end_rotation_speedup * rotation
            xtax.azimuth[] += end_rotation_speedup * rotation
        end
    end
end

const aminoacids = collect("ACDEFGHIKLMNPQRSTVWYX")
aa_to_x(aa) = findfirst(==(aa), aminoacids) / length(aminoacids)
function backbone_only!(struc)
    for chain in struc
        for atoms in chain.atoms
            deleteat!(atoms, 4:length(atoms))
        end
    end
    return struc
end

normalize_to_unit(nums) = begin
    if isempty(nums)
        Float32[]
    else
        nmin = float(minimum(nums))
        nmax = float(maximum(nums))
        denom = nmax - nmin
        denom == 0 ? fill(0.0f0, length(nums)) : [(float(x) - nmin) / denom for x in nums] .|> Float32
    end
end

function get_colors(structure, color_by::Symbol)
    if color_by == :chain
        nchains = length(structure)
        denom = max(1, nchains - 1)
        return [fill((j - 1) / denom, length(chain)) for (j, chain) in enumerate(structure)]
    elseif color_by == :aa
        return [[aa_to_x(aa) for aa in chain.sequence] for chain in structure]
    elseif color_by == :numbering
        return [normalize_to_unit(chain.numbering) for chain in structure]
    else
        throw(ArgumentError("Unsupported color_by: $(color_by)"))
    end
end

extension_whitelist(whitelist...) = f -> last(splitext(f)) in whitelist

"""
    animate_trajectory_dir(export_path, input_dirs; kwargs...)

# Keyword arguments:
- `labels`: A vector of labels for each input directory.
- `color_by`: `aa` for amino acid color, `chain` for chain color, `numbering` for numbering color.
- `aa_colormap`: A colormap for structures colored by amino acid.
- `chain_colormap`: A colormap for structures colored by chain.
- `numbering_colormap`: A colormap for structures colored by numbering.
- `rotation = 0.02`: The rotation speed.
- `framerate = 24`: The framerate of the animation.
- `end_rotation_speedup = 0.0`: Additive end rotation speed.
- `end_ribbon_seconds = 3.0`: The number of seconds to show the ribbon after the animation.
- `size`: The size of the figure.
"""
function animate_trajectory_dir(
    export_path, input_dirs;
    labels = [basename(dirpath) for dirpath in input_dirs],
    color_by = fill(:aa, length(input_dirs)), chain_colormap = :jet, aa_colormap = :turbo, numbering_colormap = :rainbow,
    rotation = 0.02, rotation_offset = 0.0, end_rotation_speedup = 0, end_ribbon_seconds = 3.0, size = (1280, 720), framerate = 24,
    theme = :black, kwargs...
)
    fig = Makie.Figure(; size, figure_padding = 1)
    theme == :black && set_theme!(theme_black())
    timestep = Observable(1)

    dir_files = Dict(dirpath => filter(extension_whitelist(".pdb", ".cif"), readdir(dirpath, join=true)) for dirpath in input_dirs)
    structures_all = Dict(dirpath => [backbone_only!(read(f, ProteinStructure)) for f in files] for (dirpath, files) in dir_files)
    structures = Dict(dirpath => @lift structures_all[dirpath][clamp($timestep, begin, end)] for (dirpath, files) in dir_files)

    dir_bounds = Dict()
    for dirpath in input_dirs
        mins = fill(typemax(Float32), 3)
        maxs = fill(typemin(Float32), 3)
        for s in structures_all[dirpath]
            for chain in s
                bb = get_backbone(chain)
                ex = extrema(bb; dims=(2, 3))
                mins = min.(mins, first.(ex))
                maxs = max.(maxs, last.(ex))
            end
        end
        dir_bounds[dirpath] = (mins[1], maxs[1], mins[2], maxs[2], mins[3], maxs[3])
    end

    axes = Any[]
    atomplot_refs = Any[]  # store plot objects to delete later
    ribbon_added = Ref(false)
    atoms_deleted = Ref(false)
    for (i, (dirpath, label)) in enumerate(zip(input_dirs, labels))
        ax = Axis3(fig[1:5, i],
            perspectiveness=0.2, protrusions = (0,0,0,0), aspect=:data, viewmode=:fit,
            width = 1000, height = 1000, tellwidth = false, tellheight = false; kwargs...)
        ax.azimuth[] = -0.2
        ax.elevation[] = 0.0
        colormap = color_by[i] == :chain ? chain_colormap : color_by[i] == :aa ? aa_colormap : numbering_colormap
        color = @lift reduce(vcat, get_colors($(structures[dirpath]), color_by[i])) |> triplicate
        ribbon_colors = @lift get_colors($(structures[dirpath]), color_by[i])
        full_label = @lift "$(basename(dirpath)): $((first ∘ splitext ∘ basename)(dir_files[dirpath][clamp($timestep, begin, end)]))"
        Label(fig[1:5, i, Top()], full_label, valign = :bottom, font = :bold, fontsize = 14, padding = (0, 0, 0, 0))
        plotref = atomplot!(ax, structures[dirpath]; color, colormap, colorrange = (0.0f0, 1.0f0))
        hidespines!(ax)
        hidedecorations!(ax)
        limits!(ax, dir_bounds[dirpath]...)
        push!(axes, ax)
        push!(atomplot_refs, (; ax, plotref, dirpath, colormap, ribbon_colors))
    end

    colgap!(fig.layout, 0)
    rowgap!(fig.layout, 0)
    maxlen = maximum(length, values(dir_files))
    end_extra_frames = Int(round(end_ribbon_seconds * framerate))
    for ax in axes
        ax.azimuth[] += rotation
    end
    record(fig, export_path, -framerate:maxlen+framerate+end_extra_frames; framerate) do t
        timestep[] = t
        for ax in axes
            ax.azimuth[] += rotation
        end
        if t >= maxlen
            for ax in axes
                ax.azimuth[] += end_rotation_speedup * rotation
            end
        end
        # Add ribbons one frame before removing atoms to reduce snapping
        if (t == maxlen + 1) && !ribbon_added[]
            for ref in atomplot_refs
                ribbon!(ref.ax, structures[ref.dirpath]; colors=ref.ribbon_colors, colormap=ref.colormap, colorrange=(0.0f0, 1.0f0))
            end
            ribbon_added[] = true
        end
        if (t == maxlen + 2) && ribbon_added[] && !atoms_deleted[]
            for ref in atomplot_refs
                delete!(ref.ax, ref.plotref)
            end
            atoms_deleted[] = true
        end
    end
end

function read_xyz_atom_line(str)
    parts = split(str)
    name = uppercase(parts[1])
    x, y, z = parse(Float64, parts[2]), parse(Float64, parts[3]), parse(Float64, parts[4])
    return Atom(name, name, [x, y, z])
end

function read_xyz(path)
    lines = readlines(path)
    atom_lines = if isnothing(tryparse(Int, lines[1]))
        filter(!isempty, lines)
    else
        natoms = parse(Int, lines[1])
        lines[3:(2+natoms)]
    end
    return read_xyz_atom_line.(atom_lines)
end

"""
    animate_molecules(export_path, molecules_cache; kwargs...)

# Keyword arguments:
- `labels`: A vector of labels for each molecule series.
- `frame_names`: A vector of vectors of frame names for each molecule series.
- `rotation = 0.02`: The rotation speed.
- `extra_seconds = 5`: The number of seconds to show the extra frames after the animation.
- `end_rotation_speedup = 0.5`: The end rotation speedup.
- `size = (1280, 720)`: The size of the figure.
- `framerate = 60`: The framerate of the animation.
- `duration = 10.0`: The duration of the animation.
- `theme = :black`: The theme of the animation.
"""
function animate_molecules(
    export_path, molecules_cache::Vector{Vector{Vector{Atom{Float64}}}};
    labels = map(string, 1:length(molecules_cache)),
    frame_names = nothing,
    rotation = 0.02, extra_seconds = 5, end_rotation_speedup = 0.5, size = (1280, 720),
    framerate = 60, duration = 10.0,
    theme = :black, kwargs...
)
    length(labels) == length(molecules_cache) || throw(ArgumentError("length(labels) must match number of molecule series"))
    if frame_names !== nothing
        length(frame_names) == length(molecules_cache) || throw(ArgumentError("length(frame_names) must match number of molecule series"))
        for i in 1:length(molecules_cache)
            length(frame_names[i]) == length(molecules_cache[i]) || throw(ArgumentError("frame_names[i] length must match molecules_cache[i] length for all i"))
        end
    end

    fig = Makie.Figure(; size, figure_padding = 1)
    theme == :black && set_theme!(theme_black())
    timestep = Observable(1)

    current_frame = [@lift molecules_cache[i][clamp($timestep, begin, end)] for i in 1:length(molecules_cache)]

    bounds = NTuple{6,Float32}[]
    for mols in molecules_cache
        mins = fill(typemax(Float32), 3)
        maxs = fill(typemin(Float32), 3)
        for molecule in mols
            coords = stack(atom_coords.(molecule))
            ex = extrema(coords; dims=2)
            mins .= min.(mins, first.(ex))
            maxs .= max.(maxs, last.(ex))
        end
        push!(bounds, (mins[1], maxs[1], mins[2], maxs[2], mins[3], maxs[3]) .* 1.6)
    end

    axes = Any[]
    for i in 1:length(molecules_cache)
        ax = Axis3(fig[1:5, i],
            perspectiveness=0.2, protrusions = (0,0,0,0), aspect=:data, viewmode=:fit,
            width = 1000, height = 1000, tellwidth = false, tellheight = false; kwargs...)
        ax.azimuth[] = -0.2
        ax.elevation[] = 0.0
        push!(axes, ax)
        if frame_names === nothing
            Label(fig[1:5, i, Top()], labels[i], valign = :bottom, font = :bold, fontsize = 14, padding = (0, 0, 0, 0))
        else
            full_label = @lift "$(labels[i]): $(frame_names[i][clamp($timestep, 1, length(frame_names[i]))])"
            Label(fig[1:5, i, Top()], full_label, valign = :bottom, font = :bold, fontsize = 14, padding = (0, 0, 0, 0))
        end
        atomplot!(ax, current_frame[i], default_size = 0.5f0, show_bonds = true, bond_width = 0.2f0)
        hidespines!(ax)
        hidedecorations!(ax)
        limits!(ax, bounds[i]...)
    end

    colgap!(fig.layout, 0)
    rowgap!(fig.layout, 0)
    maxlen = maximum(length, molecules_cache)
    end_extra_frames = Int(round(extra_seconds * framerate))
    play_frames = max(1, Int(round(framerate * duration)))
    record(fig, export_path, -framerate:play_frames+end_extra_frames; framerate) do t
        if t <= 0
            timestep[] = 1
        elseif t <= play_frames
            if play_frames <= 1
                timestep[] = 1
            else
                idx = 1 + round(Int, (t - 1) * (maxlen - 1) / (play_frames - 1))
                timestep[] = clamp(idx, 1, maxlen)
            end
        else
            timestep[] = maxlen
        end
        for ax in axes
            ax.azimuth[] += rotation
        end
        if t >= play_frames
            for ax in axes
                ax.azimuth[] += end_rotation_speedup * rotation
            end
        end
    end
end

"""
    animate_molecule_dir(export_path, input_dirs; kwargs...)

Reads all .xyz files in each input directory and animates them using `animate_molecules`.
"""
function animate_molecule_dir(
    export_path, input_dirs; labels = [basename(dirpath) for dirpath in input_dirs], kwargs...
)
    length(labels) == length(input_dirs) || throw(ArgumentError("length(labels) must match number of input directories"))
    files_per_dir = [filter(extension_whitelist(".xyz"), readdir(dirpath, join=true)) for dirpath in input_dirs]
    molecules_cache = [ [read_xyz(f) for f in files] for files in files_per_dir ]
    frame_names = [ [(first ∘ splitext ∘ basename)(f) for f in files] for files in files_per_dir ]
    animate_molecules(export_path, molecules_cache; labels, frame_names, kwargs...)
end

using ProteinChains: atom_name, atom_symbol, atom_coords

shift_atom(atom, v) = Atom(atom_name(atom), atom_symbol(atom), atom_coords(atom) .+ v)
shift_molecule(mol, v) = map(atom -> shift_atom(atom, v), mol)

function static_trajectory(input_dir;
    shift_vector=Vec3f(50.0,0,0), # offset of final molecule
    mask_color=:orange,
    size=(1920, 1080),
    alpha=0.3,
    markersize=160f0,
    kwargs...
)
    files = filter(extension_whitelist(".xyz"), readdir(input_dir, join=true))
    molecules = [read_xyz(f) for f in files]
    shift_vector_step = shift_vector / length(molecules)
    shifted_molecules = [shift_molecule(mol, i * shift_vector_step) for (i, mol) in enumerate(molecules)]
    all_atoms = reduce(vcat, shifted_molecules)
    color = [get(ATOM_COLORS, sym, mask_color) for sym in atom_symbol.(all_atoms)]
    set_theme!(theme_black())
    fig = Figure(; size)
    ax = Axis3(fig[1,1], aspect=:data)
    scatter!(ax, atom_coords.(all_atoms); color, fxaa=true, alpha, markersize, transform_marker=true, depthsorting=true)
    atomplot!(ax, shifted_molecules[end]; default_size = 0.5f0, show_bonds = true, bond_width = 0.2f0, mask_color)
    hidespines!(ax)
    hidedecorations!(ax)
    fig
end
