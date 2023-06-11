
function cumulative_digit_count(N::Integer)
    digit_sum = 0

    digits_per_number = 1
    milestone = 10
    for i in 1:N
        if i == milestone
            digits_per_number += 1
            milestone *= 10
        end
        digit_sum += digits_per_number
    end

    digit_sum
end


function find_subarrays(arr::Array{Bool,1}, L::Int, T::Int)
    subarrays::Vector{Tuple{Int,Int,Int}} = []
    N = length(arr)
    if N == 0 return subarrays end
    true_count = sum(arr[1:min(L, end)])
    for i in 2:N-L+1
        true_count += arr[i+L-1] - arr[i-1]
        if true_count >= T
            push!(subarrays, (i, i+L-1, true_count))
        end
    end
    return subarrays
end


function remove_overlapping_pairs(pairs::Vector{Tuple{Int,Int,Int}})
    N = length(pairs)
    filtered_pairs::Vector{Tuple{Int,Int,Int}} = []
    if N == 0 return filtered_pairs end
    pairs = sort(pairs, by=x -> x[1])
    current_pair = pairs[1]
    for i in 2:N
        next_pair = pairs[i]
        if next_pair[1] > current_pair[2]
            push!(filtered_pairs, current_pair)
            current_pair = next_pair
        else
            if next_pair[3] > current_pair[3]#next_pair[2]-next_pair[1] > current_pair[2]-current_pair[1]
                current_pair = next_pair
            end
        end
    end
    push!(filtered_pairs, current_pair)
    return filtered_pairs
end


"""
Longest increasing subsequence
"""
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


"""
Longest increasing subsequence elements
"""
function LISE(arr::Vector{Int})
    n = length(arr)
    dp = fill(1, n)
    parents = fill(-1, n)
    max_len, max_end = 0, -1
    for i in 1:n
        for j in 1:i
            if arr[i] > arr[j]
                if dp[j] + 1 > dp[i]
                    dp[i] = dp[j] + 1
                    parents[i] = j
                    if dp[i] > max_len
                        max_len = dp[i]
                        max_end = i
                    end
                end
            end
        end
    end
    start_index = max_end
    while start_index > 0 && parents[start_index] != -1
        start_index = parents[start_index]
    end
    result = zeros(Int, max_len)
    current_index = max_end
    for i in max_len:-1:1
        result[i] = arr[current_index]
        current_index = parents[current_index]
    end
    return result
end

export LISE

function longest_span_sequence(sorted_array::Vector{Int}, span::Int)
    if isempty(sorted_array) || span < 0
        return Vector{Int}()
    end

    start_index = 1
    end_index = 1
    longest_start = 1
    longest_end = 1

    while end_index <= length(sorted_array)
        current_span = sorted_array[end_index] - sorted_array[start_index]

        if current_span <= span
            if (end_index - start_index) >= (longest_end - longest_start)
                longest_start = start_index
                longest_end = end_index
            end
            end_index += 1
        else
            start_index += 1
        end
    end

    return sorted_array[longest_start:longest_end]
end


export longest_span_sequence


function nunique(a::SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true})
    last = first(a)
    n = 1
    for x in a
        if isless(last, x)
            n += 1
            last = x
        end
    end
    n
end

export nunique

"""
Finds the subsequence with the highest score (based on score_function),
whose first and last elements have at most a difference of window_size
"""
function max_subseq_in_range(vector::Vector{Int64}, window_size::Int64, score_function::Function=nunique)
    len = length(vector)
    current_end::Int64 = 0
    best_start::Int64 = 0
    best_end::Int64 = 0
    best_score = 0
    for current_start::Int64 in 1:len
        while current_end+1 <= len && vector[current_end+1] - vector[current_start] < window_size
            current_end += 1
        end
        current_score = score_function(view(vector, current_start:current_end))
        if current_score > best_score
            best_start = current_start
            best_end = current_end
            best_score = current_score
        end
    end
    return best_score, best_start, best_end
end

export max_subseq_in_range



function increment_dict_value!(
    key,
    delta,
    dict::Dict,
)
    if haskey(dict, key)
        dict[key] += delta
    else
        dict[key] = delta
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


@inline function constrained_LIS(arr::Vector{Int}, span::Int)
    n = length(arr)
    if n == 0
        return (0, 0, 0)
    end

    dp = ones(Int, n)
    max_lis = 1
    max_lis_end = 1

    for i in 2:n
        for j in 1:i-1
            if arr[i] > arr[j] && arr[i] - arr[j] <= span
                if dp[j] + 1 > dp[i]
                    dp[i] = dp[j] + 1
                    if dp[i] > max_lis
                        max_lis = dp[i]
                        max_lis_end = i
                    end
                end
            end
        end
    end

    max_lis_start = max_lis_end
    for i in max_lis_end-1:-1:1
        if dp[i] == dp[max_lis_start] - 1 && arr[i] < arr[max_lis_start] && arr[max_lis_end] - arr[i] <= span
            max_lis_start = i
        end
    end

    return (max_lis, arr[max_lis_start], arr[max_lis_end])
end

@inline function constrained_multi_choice_LIS(vectors::Vector{Vector{Int}}, span::Int)
    # Reverse sort each vector
    for vec in vectors
        sort!(vec, rev=true)
    end

    # Concatenate all vectors into one
    concatenated = reduce(vcat, vectors)

    lise = LISE(concatenated)

    validated_lise = longest_span_sequence(lise, span)

    # Find the longest increasing subsequence with the span constraint
    return length(validated_lise), validated_lise[1], validated_lise[end]
end
