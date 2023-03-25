

function snipe_reads(
    reference_file_path::AbstractString,
    datafiles::Vector{<:AbstractString};
    datafile_dir::AbstractString = "",
    output_dir::AbstractString = "",
    k::Int = 7,
    step::Int = 1,
    write_reads::Bool = false,
    save_data::Bool = true,
    create_plots::Bool = true,
    show_progress::Bool = true,
    t_adjustment::Int = 0,
)

    reference = Reference(reference_file_path, k)
    showinfo(reference)

    datafile_paths = joinpath.(datafile_dir, datafiles)
    datafile_list = DatafileMetadata(datafile_paths)
    
    for (i, datafile) in enumerate(datafile_list)
        showinfo(i, datafile)
    end

    # calculate read score threshold for filtering 
    expected_score_distribution = estimate_score_distribution(k, reference, datafile_list)
    threshold = max(
        first_true_after_false(x->x<1, expected_score_distribution),
        findfirst(==(maximum(expected_score_distribution)), expected_score_distribution)
    ) + t_adjustment

    # set up a config that includes run parameters
    nthreads = Threads.nthreads()
    config = Config(k, step, threshold, nthreads)
    @info "Config" k step threshold nthreads

    if show_progress println() end
    reads, score_bins = analyse_dataset(
        config,
        reference,
        datafile_list,
        show_progress,
    )
    if show_progress println() end

    if write_reads
        reads_file_basename = split(splitext(basename(datafiles[1]))[1], "_")[1]
        reads_file_path = joinpath(output_dir, reads_file_basename)*"-sniped.fasta"
        @info "Writing sniped reads to '$(reads_file_path)'"
        
        read_index_set = Set(reads.read_index)
        filter_fasta(read_index_set, datafile_paths, reads_file_path)
    end

    score_vector = getindex.(score_bins, 1)
    frequency_vector = getindex.(score_bins, 2)
    
    suffix = "k$(k)s$(step)"
    result_dir = joinpath(output_dir, replace(splitext("$(suffix)-$(now())")[1], ":" => ""))

    if (save_data || create_plots) && !ispath(result_dir)
        mkpath(result_dir) 
    end

    if save_data
        @info "Saving files to '$(result_dir)'..."

        if !ispath(result_dir) mkpath(result_dir) end

        CSV.write(joinpath(result_dir, "reads.csv"), reads)
        CSV.write(joinpath(result_dir, "score_distribution.csv"), (score=score_vector, frequency=frequency_vector))
    end

    if create_plots
        @info "Saving plots to '$(result_dir)'..."

        score_distribution_plot(
            score_vector,
            frequency_vector,
            k,
            step,
            datafile_list[1].read_length,
            expected_score_distribution,
            joinpath(result_dir, "score_distribution.svg"),
        )

        T = 8
        thresholds = [threshold + 5*binomial(t, 2) for t in 1:T]
        upper_threshold = datafile_list[1].read_length-config.k+1
        thresholds = thresholds[thresholds .< upper_threshold]
        push!(thresholds, upper_threshold)
        
        activity_plot(
            reads,
            reference.length,
            config.k,
            thresholds,
            joinpath(result_dir, "match_activity.svg"),
            true,
        )
    end

    return reads
end

export snipe_reads


function test_snipe(k = 7, step = 3, SRR_ID::AbstractString = "SRR10873757")
    snipe_reads(
        "C:/Users/anton/RSData/reference/reference.fasta",
        ["$(SRR_ID)_1.fasta", "$(SRR_ID)_2.fasta"],
        datafile_dir = "C:/Users/anton/RSData/datasets/fasta-files/$(SRR_ID)",
        output_dir = "output",
        k = k,
        step = step,
        write_reads = true,
        save_data = true,
        create_plots = true,
    )
end

export test_snipe