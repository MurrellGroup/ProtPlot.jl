export ramachandran

function ramachandran!(ax, phi_angles, psi_angles; color=:black, kwargs...)
    scatter!(ax, rad2deg.(phi_angles), rad2deg.(psi_angles), color=color)
    return nothing
end

function ramachandran!(ax, chain::Backboner.Protein.Chain; mask_missing=true, verbose=true, kwargs...)
    bonds = Backboner.ChainedBonds(chain.backbone)
    phi_angles = Backboner.Protein.phi_angles(bonds)
    psi_angles = Backboner.Protein.psi_angles(bonds)
    if mask_missing
        contiguity_mask = chain.resnums[2:end] .== chain.resnums[1:end-1] .+ 1
        if !all(contiguity_mask)
            phi_angles = phi_angles[contiguity_mask]
            psi_angles = psi_angles[contiguity_mask]
            verbose && @warn "Discarding $(count(!, contiguity_mask)) out of $(length(contiguity_mask)) points in Ramachandran plot due to missing residues."
        end
    end
    ramachandran!(ax, phi_angles, psi_angles; verbose, kwargs...)
    return nothing
end

function ramachandran!(ax, chains::AbstractVector{Backboner.Protein.Chain}; kwargs...)
    for chain in chains
        ramachandran!(ax, chain; kwargs...)
    end
    return nothing
end

function ramachandran(x;
    size=(600, 600), title="Ramachandran Plot", xlabel="Phi", ylabel="Psi", kwargs...
)
    fig = Figure(size=size)
    ax = Axis(fig[1, 1], aspect=AxisAspect(1),
        title=title, xlabel=xlabel, ylabel=ylabel,
        xticks=[-180, -90, 0, 90, 180],
        yticks=[-180, -90, 0, 90, 180],
        limits=(-180, 180, -180, 180))
    ramachandran!(ax, x; kwargs...)
    return fig
end
