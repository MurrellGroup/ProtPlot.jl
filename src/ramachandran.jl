export ramachandran

function ramachandran!(ax, phi_angles, psi_angles; color=:black, kwargs...)
    scatter!(ax, rad2deg.(phi_angles), rad2deg.(psi_angles), color=color)
    return nothing
end

function ramachandran!(ax, chain::Backboner.Protein.Chain; color=:black, kwargs...)
    bonds = Backboner.ChainedBonds(chain.backbone)
    phi_angles = Backboner.Protein.phi_angles(bonds)
    psi_angles = Backboner.Protein.psi_angles(bonds)
    ramachandran!(ax, phi_angles, psi_angles; color=color, kwargs...)
    return nothing
end

function ramachandran!(ax, chains::AbstractVector{Backboner.Protein.Chain}; kwargs...)
    for chain in chains
        ramachandran!(ax, chain; kwargs...)
    end
    return nothing
end

function ramachandran(x;
    title="Ramachandran Plot", xlabel="Phi", ylabel="Psi", kwargs...
)
    fig = Figure()
    ax = Axis(fig[1, 1],
        title=title, xlabel=xlabel, ylabel=ylabel,
        xticks=[-180, -90, 0, 90, 180],
        yticks=[-180, -90, 0, 90, 180],
        limits=(-180, 180, -180, 180))
    ramachandran!(ax, x; kwargs...)
    return fig
end
