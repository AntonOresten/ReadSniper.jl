
using BioSequences


function create_kmer_vector(ref::LongDNASeq, k::Int64)
    return LongDNASeq[ref[i:i+k-1] for i in 1:length(ref)-k+1]
end

function create_multi_match_kmer_dict(kmer_vector::Vector{LongDNASeq})
    kmer_dict = Dict{LongSequence{DNAAlphabet{4}}, Vector{Int64}}()
    for (i, kmer) in enumerate(kmer_vector)
        if haskey(kmer_dict, kmer)
            push!(kmer_dict[kmer], i)
        else
            kmer_dict[kmer] = [i]
        end
    end

    return kmer_dict
end

function get_query_index_vectors(
        query::LongDNASeq,
        kmer_dict::Dict{LongDNASeq, Vector{Int64}},
        k::Int64,
        step::Int64)
    
    index_vectors::Vector{Vector{Int64}} = [get(kmer_dict, query[i:i+k-1], Int64[]) for i in 1:step:length(query)-k+1]
end


function single_match_kmer_dict(kmer_vector::Vector{LongDNASeq})
    kmer_dict = Dict{LongDNASeq, Int64}()
    for (i, kmer) in enumerate(kmer_vector)
        if !haskey(kmer_dict, kmer)
            kmer_dict[kmer] = i
        end
    end
    return kmer_dict
end

function get_single_match_indices(
        seq::LongDNASeq,
        kmer_dict::Dict{LongDNASeq, Int64},
        k::Int64, step::Int64=1)
    index_vector::Vector{Int64} = [get(kmer_dict, seq[i:i+k-1], 0) for i in 1:step:length(seq)-k+1]
    non_zero_indices::Vector{Int64} = zeros(Int64, count(i->(i>0), index_vector))
    i = 0
    for index in index_vector
        if index > 0
            i += 1
            non_zero_indices[i] = index
        end
    end
    return non_zero_indices
end

kmer_dict = single_match_kmer_dict(create_kmer_vector(dna"ATACATACATACATACGATCGACCTACGACTC", 4))
indices = get_single_match_indices(dna"ATACAGATCGACCATCG", kmer_dict, 4)
println(kmer_dict)
println(indices)