

function snipe_reads(
    reference_file_path::AbstractString,
    datafiles::Vector{<:AbstractString};
    datafile_dir::AbstractString = "",
    output_dir::AbstractString = "",
    k::Int64 = 7,
    step::Int64 = 1,
    show_progress::Bool = true,
    save_files::Bool = true,
    create_plots::Bool = true,
)

    reference = Reference(reference_file_path, k)
    showinfo(reference)

    datafile_paths = joinpath.(datafile_dir, datafiles)
    datafile_list = DatafileMetadata(datafile_paths)
    
    for (i, datafile) in enumerate(datafile_list)
        showinfo(i, datafile)
    end

    # calculate read score threshold for filtering 
    estimated_score_distribution = spurious_score_distribution(k, reference, datafile_list)
    threshold = max(
        first_true_after_false(x->x<1, estimated_score_distribution),
        findfirst(==(maximum(estimated_score_distribution)), estimated_score_distribution)
    )

    # set up a config that includes run parameters
    nthreads = Threads.nthreads()
    config = Config(k, step, threshold, nthreads)
    @info "Config" k step threshold nthreads
    
    suffix = "k$(k)s$(step)"
    score_distr_filename = joinpath(output_dir, "score_distr-$suffix")
    activity_filename = joinpath(output_dir, "activity-$suffix")

    println()
    reads_SOA, match_frequencies = analyse_dataset(config, reference, datafile_list, show_progress)
    println()

    if save_files
        @info "Saving files to '$(output_dir)'..."

        if output_dir != ""
            mkpath(output_dir)
            @assert ispath(output_dir)
        end

        CSV.write(joinpath(output_dir, "reads-$suffix.csv"), reads_SOA)
        CSV.write(score_distr_filename*".csv", match_frequencies)
    end
    
    # reads_SOA = CSV.read("$output_dir/reads-$suffix.csv", NamedTuple)
    # match_frequencies = CSV.read("$output_dir/score_distr-$suffix.csv", NamedTuple)

    if create_plots
        @info "Saving plots to '$(output_dir)'..."

        score_distribution_plot(
            score_distr_filename*".csv",
            k,
            step,
            datafile_list[1].read_length,
            estimated_score_distribution,
            score_distr_filename*".svg",
        )

        T = 8
        thresholds = [threshold + 5*binomial(t, 2) for t in 1:T]
        activity_plot(
            reads_SOA,
            reference.length,
            config.k,
            thresholds,
            activity_filename*".svg",
            true,
        )
    end

    return reads_SOA
end


export snipe_reads