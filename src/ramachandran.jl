export ramachandran

function ramachandran!(ax, phi_angles, psi_angles; color=:black, kwargs...)
    scatter!(ax, rad2deg.(phi_angles), rad2deg.(psi_angles), color=color)
    return nothing
end

function ramachandran!(ax, chain::Backboner.Protein.Chain; verbose=false, kwargs...)
    bonds = Backboner.ChainedBonds(chain.backbone)
    valid_residues = [false; chain.resnums[1:end-2] .+ 1 .== chain.resnums[2:end-1] .&& chain.resnums[2:end-1] .+ 1 .== chain.resnums[3:end]; false]
    phi_angles = Backboner.Protein.phi_angles(bonds)[valid_residues[1:end-1]]
    psi_angles = Backboner.Protein.psi_angles(bonds)[valid_residues[2:end]]
    verbose && count(!, valid_residues) > 2 && @warn "Discarding $(count(!, valid_residues) - 2) out of $(length(valid_residues)) residues in Ramachandran plot due to missing residues."
    ramachandran!(ax, phi_angles, psi_angles; verbose, kwargs...)
    return nothing
end

function ramachandran!(ax, chains::AbstractVector{Backboner.Protein.Chain}; kwargs...)
    for chain in chains
        ramachandran!(ax, chain; kwargs...)
    end
    return nothing
end

function ramachandran(x...;
    size=(600, 600), title="Ramachandran Plot", xlabel="Phi", ylabel="Psi", kwargs...
)
    fig = Figure(size=size)
    ax = Axis(fig[1, 1], aspect=AxisAspect(1),
        title=title, xlabel=xlabel, ylabel=ylabel,
        xticks=[-180, -90, 0, 90, 180],
        yticks=[-180, -90, 0, 90, 180],
        limits=(-180, 180, -180, 180))
    ramachandran!(ax, x...; kwargs...)
    return fig
end

ramachandran(pdbfile::AbstractString; title=basename(pdbfile), kwargs...) = ramachandran(readpdb(pdbfile); title, kwargs...)