

function snipe_reads(
    reference_file_path::AbstractString,
    fasta_files::Vector{<:AbstractString},
    fasta_files_dir::AbstractString = "";
    output_dir::AbstractString = "",
    k::Int64 = 7,
    step::Int64 = 1,
)

    reference = Reference(reference_file_path, k)

    @info """Reference sequence
        length: $(reference.length)nt,
        unique k-mers: $(reference.unique_kmer_count)
        GC-content: $(round(100*reference.gc_content, digits=2))%
     """

    fasta_file_paths = fasta_files_dir * "/" .* fasta_files
    fasta_metadata_list = FASTAMetadata(fasta_file_paths)
    
    for (i, fqm) in enumerate(fasta_metadata_list)
        @info """Dataset file $i: $(fqm.path)
            size: $(round(fqm.file_size/1e9, digits=3)) GB
            reads: $(fqm.read_count)
            read length: $(fqm.read_length)
            GC-content: $(round(100*fqm.gc_content, digits=2))%
         """
    end

    # calculate threshold
    estimated_score_distribution = spurious_score_distribution(k, reference, fasta_metadata_list)
    threshold = Int(first_true_after_false(x->x<1, estimated_score_distribution))

    nthreads = Threads.nthreads()
    config = Config(k, step, threshold, nthreads)
    @info "Config" k step threshold nthreads

    if true === true
        println()
        reads_SOA, match_frequencies = analyse_dataset(config, reference, fasta_metadata_list)
        println()

        @info "Saving files to '$(output_dir)'..."
        if isdir(output_dir) mkpath(output_dir) end

        suffix = "k$(k)s$(step)"

        CSV.write("$output_dir/reads-$suffix.csv", reads_SOA)
        CSV.write("$output_dir/score_distr-$suffix.csv", match_frequencies)
    else
        suffix = "k$(k)s$(step)"
        reads_SOA = CSV.read("$output_dir/reads-$suffix.csv", NamedTuple)
        match_frequencies = CSV.read("$output_dir/score_distr-$suffix.csv", NamedTuple)
    end

    @info "Saving plots to '$(output_dir)'..."

    score_distribution_plot(
        "$output_dir/score_distr-$suffix.csv",
        fasta_metadata_list[1].read_length,
        step,
        estimated_score_distribution,
        "$output_dir/score_distr-$suffix.svg"
    )

    T = 7
    thresholds = [threshold + 5*binomial(t, 2) for t in 1:T]
    activity_plot(
        reads_SOA,
        reference.length,
        config.k,
        thresholds,
        "$output_dir/activity-$suffix.svg",
        true,
    )
end


export snipe_reads