# # Ribbon

# ## Interactive scene

using GLMakie, ProtPlot
using ProteinChains

ribbon(pdb"1ASS")

#-

ribbon(pdb"1ASS", colors=rand.(length.(pdb"1ASS")))

# ## Add to existing scene

fig = Figure()
ax = LScene(fig[1,1])

ribbon!(ax, get_backbone(pdb"1EYE"A))
ribbon!(ax, get_backbone(pdb"1EYE"A) .+ [-30, 30, 0], colormap=:blues)

fig

# ## Attributes

#=
You may customize the geometry of the ribbon by specifying the value of attributes
in the keyword arguments of your call. Here's a list of available attributes and their defaults:
- `secondary_structures = nothing` (gets assigned by an algorithm by default;
needs to be a vector of `Vector{Int}` where `1` means loop, `2` means helix, and `3` means strand)
- `colors = nothing` (gets assigned `range(0, 1, L)` for each chain by default,
mapping to `colormap`; overrides colormap if colorants are given)
- `colormap = :jet` (see the [ColorSchemes.jl catalogue](https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/);
can also be a vector of colorants)

- `coil_diameter = 0.4`
- `coil_spline_quality = 20`
- `coil_slice_quality = 20`

- `helix_width = 2.0`
- `helix_thickness = 0.5`
- `helix_spline_quality = 20`
- `helix_slice_quality = 20`

- `strand_width = 2.0`
- `strand_thickness = 0.5`
- `strand_arrow_head_length = 5.0`
- `strand_arrow_head_width = 3.5`
- `strand_spline_quality = 20`
=#

ribbon(pdb"1ASS", 
    coil_diameter=0.6, helix_width=1.5,
    strand_width=1.5, strand_arrow_head_length=3.0)

# ## Animation

#=
Ribbon plots are rendered as a set of surfaces,
but observables don't properly propagate the updates since the number of surfaces may vary,
unlike other plot types, which can render everything in a single Makie call.
For animating ribbon plots, one solution is to delete the plot from the scene and redraw it every frame.
=#

fig = Figure()
ax = LScene(fig[1,1])

p = ribbon!(ax, pdb"1ASS")

delete!(ax.scene, p)

ribbon!(ax, pdb"1EYE")

fig