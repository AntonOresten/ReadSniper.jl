
using BioSequences
using FASTX


function lis(vec::Vector{Int64})
    n = length(vec)
    lis_values = zeros(Int64, n)

    for i in 2:n
        for j in 1:i
            if vec[i] > vec[j] && lis_values[i] < lis_values[j] + 1
                lis_values[i] = lis_values[j] + 1
            end
        end
    end

    return maximum(lis_values) + 1
end

function max_window(vector::Vector{Int64}, window_size::Int64, score_function::Function=lis)
    len = length(vector)
    current_end::Int64 = 0
    best_start::Int64 = 0
    best_end::Int64 = 0
    best_score::Real = 0
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
    return (score=best_score, window_range=1:10)
end

#=function vcat2(vec_vec::Vector{Vector{Int64}})
    n = sum(length, vec_vec)
    concat_vec = zeros(Int64, n)

    i = 0
    for vec in vec_vec
        for v in vec
            i += 1
            concat_vec[i] = v
        end
    end
    return concat_vec
end=#

#=function longest_increasing_subsequence_mm(vector::Vector{Vector{Int64}})
    len = length(vector)
    largest_k = 0
    for current_start in 1:len
        k::Int64 = 1
        for current_end in current_start+1:len
            if vector[current_end] > vector[current_start]+k
                k += 1
            end
        end
        largest_k = max(largest_k, k)
    end
    return largest_k     
end=#

#=function max_window_in_vectors(
        vector::Vector{Vector{Int64}},
        window_size::Int64,
        score_function::Function=longest_increasing_subsequence)
    len = length(vector)
    current_end::Int64 = 0
    best_start::Int64 = 0
    best_end::Int64 = 0
    best_score::Real = 0
    for current_start::Int64 in 1:len
        

        if current_score > best_score
            best_start = current_start
            best_end = current_end
            best_score = current_score
        end
    end
    return best_score
end=#


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
