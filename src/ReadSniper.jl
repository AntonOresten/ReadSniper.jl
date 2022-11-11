
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
    single_match_indices,
    multi_match_kmer_dict,
    multi_match_indices,

    # reads.jl

    retrieve_reads,

    # datasets.jl

    filter_fastq


end # module