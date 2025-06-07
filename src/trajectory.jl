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
    fig = Makie.Figure(; size, backgroundcolor, figure_padding = 1)
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
