

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


function max_window(vector::Vector{Int64}, window_size::Int64, score_function::Function=lis)
    len = length(vector)
    current_end::Int64 = 0
    best_start::Int64 = 0
    best_end::Int64 = 0
    best_score::Number = 0
    for current_start::Int64 in 1:len
        while current_end+1 <= len && vector[current_end+1] - vector[current_start] < window_size
            current_end += 1
        end
        current_score = score_function(vector[current_start:current_end])
        if current_score > best_score
            best_start = current_start
            best_end = current_end
            best_score = current_score
        end
    end
    return best_score, best_start, best_end
end


function stack_index_range!(vector::Vector{Int64}, index_range::UnitRange{Int64})
    for i in index_range
        vector[i] += 1
    end
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


function increment_value(dict::Dict{Int64, Int64}, key::Int64)
    if haskey(dict, key)
        dict[key] += 1
    else
        dict[key] = 1
    end
end
