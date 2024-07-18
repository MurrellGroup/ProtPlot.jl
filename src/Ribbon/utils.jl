function expand_colors(colors::AbstractVector, N::Integer, weights::AbstractVector=ones(length(colors)))
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

# split a chain into subchains if the distance between consecutive carbonyl and nitrogen atoms is greater than a certain threshold
function get_subchain_ranges(chain_backbone::Array{T,3}, cn_distance_tolerance=2) where T <: Real
    chain_length = size(chain_backbone, 3)
    carbonyl_nitrogen_distances = sum(abs2, chain_backbone[:, 3, 1:end-1] - chain_backbone[:, 1, 2:end], dims=1) |> vec .|> sqrt
    backbone_splits = carbonyl_nitrogen_distances .> cn_distance_tolerance
    split_indices = findall(backbone_splits)
    ranges = UnitRange{Int}[]
    if isempty(split_indices)
        push!(ranges, 1:chain_length)
        return ranges
    end
    start_i = 0
    for i in split_indices
        push!(ranges, start_i+1:i)
        start_i = i
    end
    push!(ranges, start_i+1:chain_length)
    return ranges
end