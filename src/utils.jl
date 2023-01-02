

function lis(vector::Vector{<:Integer})
    n = length(vector)
    lis_values = zeros(Int64, n)
    for i in 2:n
        for j in 1:i
            if vector[i] > vector[j] && lis_values[i] < lis_values[j] + 1
                lis_values[i] = lis_values[j] + 1
            end
        end 
    end
    return maximum(lis_values) + 1
end

export lis


function increment_dict_value!(
    key,
    dict::Dict
)
    if haskey(dict, key)
        dict[key] += 1
    else
        dict[key] = 1
    end
end

export increment_dict_value!