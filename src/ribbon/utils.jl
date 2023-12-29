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