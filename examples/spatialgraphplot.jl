# # SpatialGraphPlot 

# ## Setup

using ProteinChains # reexports Backboner.Frames
using Manifolds, LinearAlgebra, Statistics

function initial_frames(frames₁::Frames)
    M = Rotations(3)
    E = Euclidean(3)
    R₁ = eachslice(frames₁.rotations, dims=3)
    t₁ = eachslice(frames₁.translations, dims=2)
    μ = mean(t₁)
    σ = std(t₁)
    R₀ = stack(rand(M, length(R₁)))
    t₀ = stack(rand(E, length(t₁))) .* σ .+ μ
    Frames(R₀, t₀)
end

function interpolate_frames(frames₀::Frames, frames₁::Frames, t::Number)
    M = Rotations(3) # rotation manifold
    E = Euclidean(3) # translation manifold
    Rₜ = stack(axes(frames₀.rotations, 3)) do i
        R₀ = frames₀.rotations[:,:,i]
        R₁ = frames₁.rotations[:,:,i]
        exp(M, R₀, t * log(M, R₀, R₁))
    end
    tₜ = stack(axes(frames₀.translations, 2)) do i
        t₀ = frames₀.translations[:,i]
        t₁ = frames₁.translations[:,i]
        exp(E, t₀, t * log(E, t₀, t₁))
    end
    return Frames(Rₜ, tₜ)
end;

# ## Animation

using GLMakie, ProtPlot, Printf
set_theme!(theme_black())

time = Observable(0.0)

fig = Figure(size=(800,600))
ax = Axis3(fig[1,1], title=(@lift "time = $(@sprintf("%.2f", $time))"),
    perspectiveness=0.2, aspect=:data, viewmode=:fit)

chain = pdb"1ASS"A
frames₁ = Frames(chain)
frames₀ = initial_frames(frames₁)
framesₜ = @lift interpolate_frames(frames₀, frames₁, $time)

k, n = 1, length(chain)
new_graph() = rand(n, n) .< k / n
graphs = [Observable(new_graph()) for _ in 1:2]
for (graph, color) in zip(graphs, [:blue, :red])
    spatialgraphplot!(ax, (@lift eachcol(($framesₜ).translations)), graph; color, alpha=0.3, linewidth=1.5)
end

p = atomplot!(ax, framesₜ;
    color=repeat(range(0, 1, n), inner=3), colormap=:jet, alpha=0.5);

#

record(fig, "graphs.mp4", enumerate(-0.2:0.01:1.5), framerate=48) do (i, t);
    if 0 < t <= 1
        time[] = t
        if i % 5 == 0
            for graph in graphs
                graph[] = new_graph()
            end
        end
    end
    if t == 1.2
        ribbon!(ax.scene, frames₁)
        delete!(ax, p)
    end
    autolimits!(ax)
end;

# ![](graphs.mp4)
