

function create_kmer_vector(seq::LongDNA{4}, k::Int)
    k -= 1
    [seq[i:i+k] for i in 1:length(seq)-k]
end

export create_kmer_vector


function kmer_index_dict(
    kmer_vector::Vector{LongDNA{4}},
    rev_comp::Bool=true,
)::Dict{LongDNA{4}, Vector{Int}}

    dict = Dict{LongDNA{4}, Vector{Int}}()
    for (index, kmer) in enumerate(kmer_vector)
        push_or_create!(kmer, index, dict)
        rev_comp ? push_or_create!(reverse_complement(kmer), -index, dict) : nothing
    end
    sort!.(values(dict), rev=true)
    dict
end

function kmer_index_dict(
    seq::LongDNA{4},
    k::Int,
    rev_comp::Bool=true,
)::Dict{LongDNA{4}, Vector{Int}}

    kmer_index_dict(create_kmer_vector(seq, k), rev_comp)
end

export kmer_index_dict


const Int_vector_empty = Int[]

function kmer_match_indices(
    seq::LongDNA{4},
    kmer_dict::Dict{LongDNA{4}, Vector{Int}},
    k::Int,
    step::Int=1,
)
    k -= 1
    Vector{Int}[get(kmer_dict, view(seq, i:i+k), Int_vector_empty) for i in 1:step:length(seq)-k]
end

export kmer_match_indices