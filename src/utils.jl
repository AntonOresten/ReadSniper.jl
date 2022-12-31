

function lis(vector::Vector{T}) where {T <: Int}
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


function create_dir(dir::String)
    if isdir(dir)
        return
    end

    parent_dir, folder = splitdir(dir)
    if !(parent_dir == "")
        create_dir(parent_dir)
    end
    mkdir(dir)
end