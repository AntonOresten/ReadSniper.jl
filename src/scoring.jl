

function compute_scores(
    config::Config,
    query_file::String,
    datafiles_paths::Vector{String},
)
    query_sequence, query_length, query_kmer_dict = fetch_query_data(query_file, config.k)
    score_distribution_dict::Dict{Int32, Int32} = Dict()

    # Read scores
    read_result_sets = parallel_scanning(config, datafiles_paths, query_kmer_dict, score_distribution_dict)
    read_result_SOA = sort!(StructArray(vcat([collect(rr_set) for rr_set in read_result_sets]...)), by=result->result.kmer_match_count, rev=true)
    # SOA: Structure-Of-Arrays

    # Read score distribution
    sorted_score_bins = sort!(collect(score_distribution_dict), by=pair->pair[1])
    scores = collect([pair[1] for pair in sorted_score_bins])
    counts = collect([pair[2] for pair in sorted_score_bins])

    max_score = maximum(scores)
    score_counts = zeros(Int32, max_score+1)
    for (score, count) in zip(scores, counts)
        # +1 because score can be 0
        score_counts[score+1] = count
    end

    return read_result_SOA, (scores=collect(Int32(0):Int32(max_score)), counts=score_counts)
end


function score_dataset(
    query_file::String,
    dataset_dir::String,
    datafiles::Vector{String};
    k::Int64=8,
    step::Int64=1,
)
    nthreads = Threads.nthreads()
    config = Config(k, step, 15, nthreads)

    datafiles_paths = (rstrip(dataset_dir, '/')*"/") .* datafiles

    read_scores_SOA, score_counts = compute_scores(config, query_file, datafiles_paths)

    return read_scores_SOA, score_counts
end