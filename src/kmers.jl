

function create_kmer_vector(seq::LongDNA, k::Int64)
    k -= 1
    [view(seq, i:i+k) for i in 1:length(seq)-k]
end

export create_kmer_vector


function kmer_index_dict(kmer_vector::Vector{LongSubSeq{DNAAlphabet{4}}})
    kmer_dict = Dict{LongSubSeq{DNAAlphabet{4}}, Vector{Int64}}()
    for (i, kmer) in enumerate(kmer_vector)
        if !haskey(kmer_dict, kmer)
            kmer_dict[kmer] = [i]
        else
            push!(kmer_dict[kmer], i)
        end
    end
    kmer_dict
end

function kmer_index_dict(seq::LongDNA{4}, k::Int64)
    kmer_index_dict(create_kmer_vector(seq, k))
end

export kmer_index_dict


const Int64_vector_empty = Int64[]

function kmer_match_indices(
    seq::LongDNA{4},
    kmer_dict::Dict{LongSubSeq{DNAAlphabet{4}}, Vector{Int64}},
    k::Int64,
    step::Int64=1,
)
    k -= 1
    filter(!isempty, Vector{Int64}[get(kmer_dict, view(seq, i:i+k), Int64_vector_empty) for i in 1:step:length(seq)-k]) 
end

export kmer_match_indices