
using BioSequences


function create_kmer_vector(seq::LongDNASeq, k::Int64)
    return LongDNASeq[seq[i:i+k-1] for i in 1:length(seq)-k+1]
end

function create_kmer_dict(kmer_vector::Vector{LongDNASeq})
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

