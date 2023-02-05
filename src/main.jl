

function snipe_reads(
    reference_file_path::AbstractString,
    datafiles::Vector{<:AbstractString};
    datafile_dir::AbstractString = "",
    output_dir::AbstractString = "",
    k::Int = 7,
    step::Int = 1,
    t_adjustment::Int = 0,
    show_progress::Bool = true,
    save_data::Bool = true,
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
    ) + t_adjustment

    # set up a config that includes run parameters
    nthreads = Threads.nthreads()
    config = Config(k, step, threshold, nthreads)
    @info "Config" k step threshold nthreads

    println()
    reads_SOA, match_frequencies = analyse_dataset(config, reference, datafile_list, show_progress)
    println()
    
    suffix = "k$(k)s$(step)"
    output_folder = joinpath(output_dir, replace(splitext("$(suffix)-$(now())")[1], ":" => ""))

    if save_data
        @info "Saving files to '$(output_folder)'..."

        if !ispath(output_folder) mkpath(output_folder) end

        CSV.write(joinpath(output_folder, "reads.csv"), reads_SOA)
        CSV.write(joinpath(output_folder, "match_activity.csv"), match_frequencies)
    end
    
    # reads_SOA = CSV.read("$output_folder/reads-$suffix.csv", NamedTuple)
    # match_frequencies = CSV.read("$output_folder/score_distr-$suffix.csv", NamedTuple)

    if create_plots
        @info "Saving plots to '$(output_folder)'..."
        
        if !ispath(output_folder) mkpath(output_folder) end

        score_distribution_plot(
            match_frequencies.matches,
            match_frequencies.frequency,
            k,
            step,
            datafile_list[1].read_length,
            estimated_score_distribution,
            joinpath(output_folder, "score_distribution.svg"),
        )

        T = 8
        thresholds = [threshold + 5*binomial(t, 2) for t in 1:T]
        upper_threshold = datafile_list[1].read_length-config.k+1
        thresholds = thresholds[thresholds .< upper_threshold]
        push!(thresholds, upper_threshold)
        
        activity_plot(
            reads_SOA,
            reference.length,
            config.k,
            thresholds,
            joinpath(output_folder, "match_activity.svg"),
            true,
        )
    end

    return reads_SOA
end

export snipe_reads


function test_snipe()
    snipe_reads(
        "C:/Users/anton/RSData/reference/dummy.fasta",
        ["SRR10873757_1.fasta", "SRR10873757_2.fasta"],
        datafile_dir = "C:/Users/anton/RSData/datasets/fasta-files/SRR10873757",
        output_dir = "output",
        k = 8,
        step = 3,
        save_data = true,
        create_plots = true,
    )
end

export test_snipe