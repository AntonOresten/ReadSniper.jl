
function size_score(sorted_subvector::Vector{Int64})
    return length(sorted_subvector)
end

function max_window(sorted_vector::Vector{Int64}, window_size::Int64, score_function::Function=size_score)
    len = length(sorted_vector)
    current_end::Int64 = 0
    best_start::Int64 = 0
    best_end::Int64 = 0
    best_score::Real = 0
    for current_start::Int64 in 1:len
        while current_end+1 <= len && sorted_vector[current_end+1] - sorted_vector[current_start] < window_size
            current_end += 1
        end
        current_score = score_function(sorted_vector[current_start:current_end])
        if current_score > best_score
            best_start = current_start
            best_end = current_end
            best_score = current_score
        end
        #=println("$current_start, $current_end");println("$(sorted_vector[current_start]), $(sorted_vector[current_end])");println("$current_score\n")=#
    end
    return best_score
end

@show max_window([1,2,5,6,9,10,11,15,16,17,20,25], 6)
