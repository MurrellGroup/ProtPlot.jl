export Ramachandran, ramachandran, ramachandran!

@recipe(Ramachandran, phi_angles, psi_angles) do scene
    Attributes(
        color = :black,
    )
end

function Makie.plot!(ramachandran::Ramachandran{<:Tuple{<:AbstractVector{<:Real}, <:AbstractVector{<:Real}}})
    phi_angles = ramachandran[1][]
    psi_angles = ramachandran[2][]
    scatter!(ramachandran, rad2deg.(phi_angles), rad2deg.(psi_angles), color=ramachandran.color[])
    return ramachandran
end

function Makie.convert_arguments(::Type{<:Ramachandran}, chain::Protein.Chain)
    bonds = ChainedBonds(chain.backbone)
    valid_residues = [false; chain.resnums[1:end-2] .+ 1 .== chain.resnums[2:end-1] .&& chain.resnums[2:end-1] .+ 1 .== chain.resnums[3:end]; false]
    phi_angles = Protein.phi_angles(bonds)[valid_residues[1:end-1]]
    psi_angles = Protein.psi_angles(bonds)[valid_residues[2:end]]
    count(!, valid_residues) > 2 && @warn "Discarding $(count(!, valid_residues) - 2) out of $(length(valid_residues)) residues in Ramachandran plot due to missing residues."
    return (phi_angles, psi_angles)
end

function Makie.convert_arguments(::Type{<:Ramachandran}, chains::Vector{Protein.Chain})
    phi_angles = Real[]
    psi_angles = Real[]
    for chain in chains
        phi_angles_chain, psi_angles_chain = Makie.convert_arguments(Ramachandran, chain)
        append!(phi_angles, phi_angles_chain)
        append!(psi_angles, psi_angles_chain)
    end
    return (phi_angles, psi_angles)
end

Makie.convert_arguments(T::Type{<:Ramachandran}, pdbfile::AbstractString) = Makie.convert_arguments(T, readpdb(pdbfile))