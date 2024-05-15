export animate_attention

function render_rotation_frame!(container,
    chain::Protein.Chain, frames::Frames, i::Int;
    plot_list=nothing,
)
    frame = frames[i]
    m = collect(frame.rotation)
    len = 2.0
    p = arrows!(container,
        fill(Point3f(chain.backbone[:, 3i-1]...), 3),
        [len * Vec3f(col...) for col in eachcol(m)],
        colormap=:seaborn_dark6,
        color=[5/12, 3/12, 1/12],
        arrowsize = Vec3f(0.6, 0.6, 1.0),
        linewidth=0.4,
        fxaa=true)
    !isnothing(plot_list) && push!(plot_list, p)
end

"""
    animate_attention(chain::Protein.Chain, attention::AbstractArray{<:Real, 3}; kwargs...)
"""
function animate_attention(
    chain::Protein.Chain, attention::AbstractArray{<:Real, 3};
    azimuth_start = 1, azimuth_end = -6, output_file::String = "attention.mp4",
    ribbon_colorscheme = :glasgow,
    attention_colorscheme = :hsv,
    end_padding = 3, grow_limits = false, from_centroid = false,
    frames_per_residue::Int = 10, framerate::Int = 30, show_rotation_frame = false,
    kwargs...
)
    frames = Frames(chain.backbone, Protein.STANDARD_TRIANGLE_ANGSTROM)
    if from_centroid
        points = frames.locations
    else
        points = Protein.alphacarbon_coords(chain)
    end

    attention = PointAttention(points, attention)

    fig = Figure();
    ax = Axis3(fig[1, 1], protrusions=(20, 20, 10, 10), perspectiveness=0.2, aspect=:data, viewmode=:fit);
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
    plot_list = Vector{AbstractPlot}()

    n = length(chain)
    k = frames_per_residue
    
    frame_indices = 2:1/k:n+end_padding*framerate/k
    azimuth(t) = (t / (last(frame_indices) - first(frame_indices))) * (azimuth_end - azimuth_start) + azimuth_start
    record(fig, output_file, frame_indices, framerate=framerate) do x
        ax.azimuth[] = azimuth(x)
        if x % 1 == 0 && x <= n
            i = round(Int, x)

            for plot in plot_list
                delete!(ax.scene, plot)
            end
            empty!(plot_list)

            subchain = Protein.Chain(@view(chain.backbone[1:3i]), ssvector=chain.ssvector[1:i])
            r = ribbon!(ax, [subchain], colors=[LinRange(0, i/n, i)], colormap=ribbon_colorscheme)
            push!(plot_list, r)

            show_rotation_frame && render_rotation_frame!(ax, chain, frames, i; plot_list)

            try
                H = size(attention.intensities, 1)
                Attention.draw_attention_slice!(ax, i, attention, threshold=0.01, colors=range(0, !isone(H), H), colormap=attention_colorscheme; plot_list, kwargs...)
            catch e
                println("Got $e while attempting to render attention for residue $i.")
            end
        end
    end
end