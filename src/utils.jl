

function find_subarrays(arr::Vector{Bool}, L::Int, T::Int)
    subarrays::Vector{Tuple{Int,Int}} = []
    n = length(arr)
    true_count = sum(arr[1:min(L, end)])
    for i in L+1:n
        true_count += arr[i] - arr[i-L]
        if true_count >= T
            push!(subarrays, (i-L+1,i))
        end
    end
    return subarrays
end

export find_subarrays


function remove_overlapping_pairs(pairs::Vector{Tuple{Int,Int}})
    filtered_pairs::Vector{Tuple{Int,Int}} = []
    N = length(pairs)
    if N == 0 return filtered_pairs end
    pairs = sort(pairs, by=x -> x[1])
    current_pair = pairs[1]
    for i in 2:N
        next_pair = pairs[i]
        if next_pair[1] > current_pair[2]
            push!(filtered_pairs, current_pair)
            current_pair = next_pair
        else
            if next_pair[2]-next_pair[1] > current_pair[2]-current_pair[1]
                current_pair = next_pair
            end
        end
    end
    push!(filtered_pairs, current_pair)
    return filtered_pairs
end

export remove_overlapping_pairs


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


function highest_score_subseq_in_window(vector::Vector{Int64}, window_size::Int64, score_function::Function=length)
    len = length(vector)
    current_end::Int64 = 0
    best_start::Int64 = 0
    best_end::Int64 = 0
    best_score = 0
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

export highest_score_subseq_in_window



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


function push_or_create!(
    key::KeyType,
    value::ValueType,
    dict::Dict{KeyType, Vector{ValueType}},
) where (KeyType) where (ValueType)
    if haskey(dict, key)
        push!(dict[key], value)
    else
        dict[key] = [value]
    end
end

export push_or_create!


function find_nth_true_region(n::Integer, func::Function, vector::Vector{<:Number})
    bools::Vector{Bool} = map(func, vector)
    k = bools[1] + 0
    if k == n
        return 1
    end
    for i in 1:length(bools)-1
        if !(!bools[i] && bools[i+1])
            continue
        end
        if (k += 1) == n
            return i+1
        end
    end
    return -1
end

function first_true_after_false(func::Function, vector::Vector{<:Number})
    if func(vector[1])
        return find_nth_true_region(2, func, vector)
    else
        return find_nth_true_region(1, func, vector)
    end
end

ge(t::Real, vector::Vector{<:Real}) = map(x->x>=t, vector)