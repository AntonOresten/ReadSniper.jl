
module ReadSniper

include("include_all.jl")

export

# utils.jl

lis,
max_window,
filter_out_empty_vectors,    

# kmers.jl

create_kmer_vector,
single_match_kmer_dict,
get_single_match_indices,
multi_match_kmer_dict,
get_multi_match_indices,

# reads

retrieve_reads_sm,
retrieve_reads_mm

end # module