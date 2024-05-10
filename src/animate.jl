module Animate

export animate_attention

using ..ProtPlot

using GLMakie
using ColorSchemes

function render_rotation_frame!(container,
    chain::Backboner.Protein.Chain, frames::Backboner.Frames, i::Int,
    rotation_frame_lightness::Float64; plots=nothing,
)
    frame = frames[i]
    m = collect(frame.rotation)
    len = 2.0
    p = arrows!(container,
        fill(Point3f(chain.backbone[:, 3i]...), 3),
        [len * Vec3f(col...) for col in eachcol(m)],
        color=[colorant"red", colorant"lime", colorant"blue"] .* rotation_frame_lightness,
        arrowsize = Vec3f(0.6, 0.6, 1.0),
        linewidth=0.4,
        fxaa=true)
    !isnothing(plots) && push!(plots, p)
end

"""
    animate_attention(chain::Backboner.Protein.Chain, attention::AbstractArray{<:Real, 3}; kwargs...)
"""
function animate_attention(
    chain::Backboner.Protein.Chain, attention::AbstractArray{<:Real, 3};
    azimuth_start = 1, azimuth_end = -6, output_file::String = "attention.mp4",
    ribbon_colorscheme = ColorSchemes.glasgow,
    attention_colorscheme = ColorSchemes.hsv,
    end_padding = 3, grow_limits = false, from_centroid = false,
    frames_per_residue::Int = 10, framerate::Int = 30, show_rotation_frame = false,
    rotation_frame_lightness = 0.5, kwargs...
)
    frames = Backboner.Frames(chain.backbone, Backboner.Protein.STANDARD_TRIANGLE_ANGSTROM)
    if from_centroid
        points = frames.locations
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
    record(fig, output_file, frame_indices, framerate=framerate) do x
        ax.azimuth[] = azimuth(x)
        if x % 1 == 0 && x <= n
            i = round(Int, x)

            for plot in plots
                delete!(ax.scene, plot)
            end
            empty!(plots)

            show_rotation_frame && render_rotation_frame!(ax, chain, frames, i, rotation_frame_lightness, plots=plots)

            subchain = Backboner.Protein.Chain(@view(chain.backbone[1:3i]), ssvector=chain.ssvector[1:i])
            ribbon!(ax, [subchain], colorscheme=ribbon_colorscheme, color_vectors=[range(0, i/length(chain), i)], plots=plots)

            try
                H = size(attention.intensities, 1)
                Attention.draw_attention_slice!(ax, i, attention, threshold=0.01, colors=attention_colorscheme[isone(H) ? [0.0] : range(0, 1, H)]; plots, kwargs...)
            catch e
                println("Got $e while attempting to render attention for residue $i.")
            end
        end
    end
end

end