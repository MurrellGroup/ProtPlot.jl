function remove_singleton_strands!(chain::Chain)
    ssvector = chain.ssvector
    for i in 2:length(ssvector)-1
        if ssvector[i-1] != Strand && ssvector[i] == Strand && ssvector[i+1] != Strand
            ssvector[i] = Loop
        end
    end
end

# expand vector to have length N while keeping the same discrete color gradient 
function expand_colors(colors::Vector, N::Integer)
    L = length(colors)
    color_matrix = reshape(repeat(colors, inner=ceil(Int, N/L))[1:N], :, 1)
    return color_matrix
end

function chain_color_vector(chain::Chain, colorscheme::ColorScheme)
    #return repeat([colorant"red", colorant"blue"], length(chain) ÷ 2 + 1)[1:length(chain)]
    return colorscheme[LinRange(0, 1, length(chain))]
end