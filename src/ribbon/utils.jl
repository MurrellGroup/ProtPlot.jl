function smooth_color_vector(colorscheme, N::Integer)
    return colorscheme[LinRange(0, 1, N)]
end

function expand_colors(colors::Vector, N::Integer)
    L = length(colors)
    repeats = ceil(Int, N / L)
    total_elements = repeats * L
    extra_elements = total_elements - N

    result = Vector{eltype(colors)}()
    for i in 1:L
        current_repeats = repeats - (i <= extra_elements ? 1 : 0)
        append!(result, fill(colors[i], current_repeats))
    end

    return reshape(result[1:N], :, 1)
end

# split a chain into subchains if there is a gap in the residue numbering or if the distance between consecutive carbonyl and nitrogen atoms is greater than a certain threshold
function Base.split(chain::Protein.Chain; resnums=true, cn_distance_tolerance=2)
    resnum_splits = resnums ? (diff(chain.resnums) .!= 1) : falses(length(chain.resnums)-1)
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