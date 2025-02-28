# # Ramachandran

using GLMakie, ProtPlot

set_theme!(theme_black())

fig = Figure()
ax = Axis(fig[1,1], title="Ramachandran plot",
    xlabel="Phi", ylabel="Psi",
    yticks=(-180:90:180, ["-180°", "-90°", "0°", "90°", "180°"]),
    xticks=(-180:90:180, ["-180°", "-90°", "0°", "90°", "180°"]),
    limits=((-180, 180), (-180, 180)))

ramachandran!(ax, pdb"1ASS", color=:white)

fig
