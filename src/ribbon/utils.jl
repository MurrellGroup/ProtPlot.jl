function smooth_color_vector(colorscheme::ColorScheme, N::Integer)
    return colorscheme[LinRange(0, 1, N)]
end

# expand vector to have length N while keeping the same discrete color scheme
function expand_colors(colors::Vector, N::Integer)
    L = length(colors)
    repeats = ceil(Int, N / L)
    total_elements = repeats * L
    extra_elements = total_elements - N

    # Distribute the extra elements evenly across the vector
    result = Vector{eltype(colors)}()
    for i in 1:L
        # Calculate the number of repeats for this element
        current_repeats = repeats - (i <= extra_elements ? 1 : 0)
        append!(result, fill(colors[i], current_repeats))
    end

    return reshape(result[1:N], :, 1)
end

function Base.getindex(chain::Protein.Chain, r::UnitRange{<:Integer})
    Protein.Chain(chain.id, chain.backbone[3*r.start-2:3r.stop], modelnum=chain.modelnum, resnums=chain.resnums[r], aavector=chain.aavector[r], ssvector=chain.ssvector[r])
end

# split a chain into subchains if there is a gap in the residue numbering
function split_by_resnum(chain::Protein.Chain, cn_distance_tolerance = 2)
    resnum_splits = diff(chain.resnums) .!= 1
    backbone_splits = Protein.carbonyl_nitrogen_distances(chain) .> cn_distance_tolerance
    split_indices = findall(resnum_splits .| backbone_splits)
    if isempty(split_indices)
        return [chain]
    end
    ranges = UnitRange{Int}[]
    start_idx = 0
    for idx in split_indices
        push!(ranges, start_idx+1:idx)
        start_idx = idx
    end
    push!(ranges, start_idx+1:length(chain))
    return ranges
end