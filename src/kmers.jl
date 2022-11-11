
function create_kmer_vector(seq::LongDNASeq, k::Int64)
    return [kmer for (_, kmer) in each(DNAMer{k}, seq)]
end

function single_match_kmer_dict(kmer_vector::Vector{DNAMer{k}}) where {k}
    kmer_dict = Dict{DNAMer{k}, Int64}()
    for (i, kmer) in enumerate(kmer_vector)
        if !haskey(kmer_dict, kmer)
            kmer_dict[kmer] = i
        end
    end
    return kmer_dict
end

function single_match_kmer_dict(seq::LongDNASeq, k::Int64)
    kmer_dict = Dict{DNAMer{k}, Int64}()
    for (i, kmer) in each(DNAMer{k}, seq)
        if !haskey(kmer_dict, kmer)
            kmer_dict[kmer] = i
        end
    end
    return kmer_dict
end

function single_match_indices(
        seq::LongDNASeq,
        kmer_dict::Dict{DNAMer{k}, Int64},
        step::Int64=1) where {k}
    index_vector::Vector{Int64} = [get(kmer_dict, kmer, 0) for (_, kmer) in each(DNAMer{k}, seq, step)]
    non_zero_indices = index_vector[index_vector .!= 0]
    return non_zero_indices
    #=non_zero_indices::Vector{Int64} = zeros(Int64, count(i->(i>0), index_vector))
    i = 0
    for index in index_vector
        if index > 0
            i += 1
            non_zero_indices[i] = index
        end
    end
    return non_zero_indices=#
end



function multi_match_kmer_dict(kmer_vector::Vector{DNAMer{k}}) where {k}
    kmer_dict = Dict{DNAMer{k}, Vector{Int64}}()
    for (i, kmer) in enumerate(kmer_vector)
        if !haskey(kmer_dict, kmer)
            kmer_dict[kmer] = [i]
        else
            push!(kmer_dict[kmer], i)
        end
    end
    return kmer_dict
end

function multi_match_kmer_dict(seq::LongDNASeq, k::Int64)
    kmer_dict = Dict{DNAMer{k}, Vector{Int64}}()
    for (i, kmer) in each(DNAMer{k}, seq)
        if !haskey(kmer_dict, kmer)
            kmer_dict[kmer] = [i]
        else
            push!(kmer_dict[kmer], i)
        end
    end
    return kmer_dict
end

function multi_match_indices(
        seq::LongDNASeq,
        kmer_dict::Dict{DNAMer{k}, Int64},
        step::Int64=1) where {k}
    index_vectors::Vector{Vector{Int64}} = [get(kmer_dict, kmer, Int64[]) for (_, kmer) in each(DNAMer{k}, seq, step)]
    non_zero_index_vectors::Vector{Vector{Int64}} = fill(Int64[], count(v->!isempty(v), index_vectors))
    i = 0
    for iv in index_vectors
        if !isempty(iv)
            i += 1
            non_zero_index_vectors[i] = iv
        end
    end
    return non_zero_index_vectors
end

#=kmer_dict = single_match_kmer_dict(create_kmer_vector(dna"ATACATACATACATACGATCGACCTACGACTC", 4))
indices = single_match_indices(dna"ATACAGATCGACCATCG", kmer_dict, 4)
println(kmer_dict)
println(indices)=#

