module Animate

export animate_attention

using ..ProtPlot

using GLMakie
using ColorSchemes

"""
    animate_attention(chain::Backboner.Protein.Chain, attention::AbstractArray{<:Real, 3}; kwargs...)
"""
function animate_attention(
    chain::Backboner.Protein.Chain, attention::AbstractArray{<:Real, 3};
    azimuth_start = 1, azimuth_end = -6, output_file::String = "attention.mp4",
    ribbon_colorscheme = ColorSchemes.glasgow,
    attention_colorscheme = ColorSchemes.hsv,
    end_padding = 3, grow_limits = false, from_centroid = true, frames_per_residue::Int = 10, framerate::Int = 30
)
    
    if from_centroid
        points = Backboner.Frames(chain.backbone,Backboner.Protein.STANDARD_TRIANGLE_ANGSTROM).locations
    else
        points = Backboner.Protein.carbonyl_coords(chain)
    end

    attention = PointAttention(points, attention)

    fig = Figure();
    ax = Axis3(fig[1, 1], protrusions=(20, 20, 10, 10), perspectiveness=0.2, aspect=:data);
    hidespines!(ax)
    hidedecorations!(ax)

    if !grow_limits
        plot_limits = extrema(points, dims=2)
        xlims!(ax, plot_limits[1]) 
        ylims!(ax, plot_limits[2]) 
        zlims!(ax, plot_limits[3]) 
    end

    ax.azimuth[] = azimuth_start

    # Aggregate all the plots from each frame, such that they can be deleted.
    # Ribbon plots are actually made up of multiple plots, so each subplot gets added to the list.
    # A possible optimization is only deleting changed segments of the ribbon plot (e.g. last segment, and coils with emerging beta sheets)
    plots = Vector{AbstractPlot}()

    n = length(chain)
    k = frames_per_residue
    
    frame_indices = 2:1/k:n+end_padding*framerate/k
    azimuth(t) = (t / (last(frame_indices) - first(frame_indices))) * (azimuth_end - azimuth_start) + azimuth_start
    record(fig, output_file, frame_indices, framerate=framerate) do i
        ax.azimuth[] = azimuth(i)
        if i % 1 == 0 && i <= n
            i = round(Int, i)

            for plot in plots
                delete!(ax.scene, plot)
            end
            empty!(plots)

            subchain = Backboner.Protein.Chain(@view(chain.backbone[1:3i]), ssvector=chain.ssvector[1:i])
            ribbon!(ax, [subchain], colorscheme=ribbon_colorscheme, color_vectors=[range(0, i/length(chain), i)], plots=plots)

            try
                H = size(attention.intensities, 1)
                Attention.draw_attention_slice!(ax, i, attention, threshold=0.01, colors=attention_colorscheme[isone(H) ? [0.0] : range(0, 1, H)], plots=plots)
            catch e
                println("Got $e while attempting to render attention for residue $i.")
            end
        end
    end
end

end