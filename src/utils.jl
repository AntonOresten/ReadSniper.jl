

function lis(vec::Vector{Int64})
    n = length(vec)
    lis_values = zeros(Int64, n)
    for i in 2:n
        for j in 1:i
            if vec[i] > vec[j] && lis_values[i] < lis_values[j] + 1
                lis_values[i] = lis_values[j] + 1
    end end end
    return maximum(lis_values) + 1
end


function filter_out_empty_vectors(vectors::Vector{Vector{Int64}})

    # Make space for non-empty vectors
    number_of_non_empty::Int64 = count(v->!isempty(v), vectors)
    non_empty_vectors::Vector{Vector{Int64}} = fill(Int64[], number_of_non_empty)
    
    i = 0
    for v in vectors
        if !isempty(v)
            i += 1
            non_empty_vectors[i] = v
        end
    end

    return non_empty_vectors
end


function increment_dict_value!(key::Int64, dict::Dict{Int64, Int64})
    if haskey(dict, key)
        dict[key] += 1
    else
        dict[key] = 1
    end
end
