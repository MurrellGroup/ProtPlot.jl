@recipe(Ribbon, chains) do scene
    Attributes(
        secondary_structures = nothing,
        colors = nothing,
        colormap = :jet,

        coil_diameter = 0.4,
        coil_spline_quality = 20,
        coil_slice_quality = 20,

        helix_width = 2.0,
        helix_thickness = 0.5,
        helix_spline_quality = 20,
        helix_slice_quality = 20,

        strand_width = 2.0,
        strand_thickness = 0.5,
        strand_spline_quality = 20,
        strand_arrow_head_length = 5,
        strand_arrow_head_width = 3.5,
    )
end

include("utils.jl")
include("secondary.jl")
include("render.jl")

# TODO: observe chains and re-render when they change
function Makie.plot!(ribbon::Ribbon{<:Tuple{<:AbstractVector{<:AbstractArray{T,3}}}}) where T<:Real
    chain_backbones = convert.(Array{T,3}, collect(ribbon[1][]))
    isnothing(ribbon.secondary_structures[]) && (ribbon.secondary_structures[] = _assign_secondary_structure(chain_backbones))
    isnothing(ribbon.colors[]) && (ribbon.colors[] = [range(0, !isone(L), L) for L in size.(chain_backbones, 3)])
    render!(ribbon, chain_backbones)
    return ribbon
end

Makie.convert_arguments(::Type{<:Ribbon}, chain_backbone::AbstractArray{T,3}) where T<:Real = ([coords],)
Makie.convert_arguments(R::Type{<:Ribbon}) = Makie.convert_arguments(R, Array{Float64,3}(undef, 3, 3, 0))

Makie.convert_arguments(R::Type{<:Ribbon}, chains::AbstractVector{<:ProteinChain}) = Makie.convert_arguments(R, map(chain -> chain.backbone, chains))
Makie.convert_arguments(R::Type{<:Ribbon}, chain::ProteinChain) = Makie.convert_arguments(R, [chain])

Makie.convert_arguments(R::Type{<:Ribbon}, path::AbstractString, args...) = Makie.convert_arguments(R, readchains(path, args...))

"""
    ribbon_scene(chains::AbstractVector{Protein.Chain}; backgroundcolor=:black, camcontrols=(;), kwargs...)

Render a protein as a ribbon diagram. The display will be automatically centered on the rendered ribbon,
unless the user supplies `camcontrols` (see Makie's camera documentation for details).
"""
function ribbon_scene(args...; backgroundcolor=:black, camcontrols=(;), kwargs...)
    scene = Scene(backgroundcolor=backgroundcolor)
    cam3d!(scene; camcontrols...)
    ribbon!(scene, args...; kwargs...)
    isempty(camcontrols) && center!(scene)
    return scene
end
