function expand_colors(
    colors::AbstractVector, N::Integer, weights::AbstractVector=ones(length(colors))
)
    L = length(colors)
    weight_sum = sum(weights)

    L != length(weights) && throw(ArgumentError("The lengths of colors and weights must be the same."))
    any(<(0), weights) && throw(ArgumentError("All weights must be non-negative."))
    iszero(weight_sum) && throw(ArgumentError("Sum of weights must not be zero."))

    normalized_weights = weights / weight_sum

    result = similar(colors, 0)
    for i in 1:L
        current_repeats = floor(Int, N * normalized_weights[i])
        append!(result, fill(colors[i], current_repeats))
    end

    while length(result) > N
        pop!(result)
    end
    while length(result) < N
        append!(result, colors[end])
    end

    return reshape(result, :, 1)
end


# split a chain into subchains if there is a gap in the residue numbering or if the distance between consecutive carbonyl and nitrogen atoms is greater than a certain threshold
function get_subchain_ranges(chain::Protein.Chain; resnums=true, cn_distance_tolerance=2)
    resnum_splits = resnums ? (diff(chain.resnums) .!= 1) : falses(length(chain.resnums)-1)
    backbone_splits = Protein.carbonyl_nitrogen_distances(chain) .> cn_distance_tolerance
    split_indices = findall(resnum_splits .| backbone_splits)
    ranges = UnitRange{Int}[]
    if isempty(split_indices)
        push!(ranges, 1:length(chain))
        return ranges
    end
    start_i = 0
    for i in split_indices
        push!(ranges, start_i+1:i)
        start_i = i
    end
    push!(ranges, start_i+1:length(chain))
    return ranges
end