module ReadSniper
    include("include_all.jl")

    export

    # utils.jl

    size_score,
    maximum_increasing_subsequence,
    max_window,    

    # kmers.jl

    create_kmer_vector,
    single_match_kmer_dict,
    get_single_match_indices

end