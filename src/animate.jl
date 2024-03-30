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
    azimuth_start = 1, azimuth_end = 0, output_file::String = "attention.mp4",
    ribbon_colorscheme = ColorSchemes.jet,
    attention_colorscheme = ColorSchemes.hawaii,
    end_padding = 3,
)
    points = Backboner.Protein.alphacarbon_coords(chain)
    attention = PointAttention(points, attention)

    fig = Figure();
    ax = Axis3(fig[1, 1], protrusions=(20, 20, 10, 10), perspectiveness=0.2, aspect=:data);
    hidespines!(ax)
    #hidedecorations!(ax)

    ax.azimuth[] = azimuth_start

    # Aggregate all the plots from each frame, such that they can be deleted.
    # Ribbon plots are actually made up of multiple plots, so each subplot gets added to the list.
    # A possible optimization is only deleting the last segment of the ribbon plot.
    plots = Vector{AbstractPlot}()

    n = length(chain)
    k = 5 # 5 frames per residue
    framerate = 30
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
                Attention.draw_attention_slice!(ax, i, attention, threshold=0.01, colors=attention_colorscheme[range(0, 1, size(attention.intensities, 1))], plots=plots)
            catch e
                println("Got $e while attempting to render attention for residue $i.")
            end
        end
    end
end

end