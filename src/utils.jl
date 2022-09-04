
using BioSequences
using FASTX

function size_score(sorted_subvector::Vector{Int64})
    return length(sorted_subvector)
end

function maximum_increasing_subsequence(vector::Vector{Int64})
    len = length(vector)
    largest_k = 0
    for current_start in 1:len
        k::Int64 = 1
        for current_end in current_start+1:len
            if vector[current_end] >= vector[current_start]+k
                k += 1
            end
        end
        largest_k = max(largest_k, k)
    end
    return largest_k     
end

function max_window(vector::Vector{Int64}, window_size::Int64, score_function::Function=maximum_increasing_subsequence)
    len = length(vector)
    current_end::Int64 = 0
    best_start::Int64 = 0
    best_end::Int64 = 0
    best_score::Real = 0
    for current_start::Int64 in 1:len
        while current_end+1 <= len && vector[current_end+1] - vector[current_start] < window_size
            current_end += 1
        end
        current_score = maximum_increasing_subsequence(vector[current_start:current_end])
        if current_score > best_score
            best_start = current_start
            best_end = current_end
            best_score = current_score
        end
        #=println("$current_start, $current_end");println("$(sorted_vector[current_start]), $(sorted_vector[current_end])");println("$current_score\n")=#
    end
    println(best_start, best_end)
    return best_score
end

@show maximum_increasing_subsequence([1,-1,2,3,-5,-4,-5,4,5,6,7])

#@show max_window([1,2,5,6,9,10,11,15,16,17,20,25], 6)


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
