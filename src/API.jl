

function match_frequency_plot(
    file::AbstractString,
    xlim_upper::Integer = 250,
    step::Integer = 1,
    outfile::AbstractString = "match_frequencies.svg",
)
    match_frequencies = CSV.read(file, NamedTuple)
    match_frequency_plot(match_frequencies.matches, match_frequencies.frequency, xlim_upper, step, outfile)
end

export match_frequency_plot


function snipe_reads(
    ref_genome_file_path::AbstractString,
    fastq_file_dir::AbstractString,
    fastq_files::Vector{<:AbstractString};
    output_dir::AbstractString="",
    k::Int64=7,
    step::Int64=1,
)
    nthreads = Threads.nthreads()
    config = Config(k, step, 30, nthreads)

    fastq_file_paths = (rstrip(fastq_file_dir, '/')*"/") .* fastq_files

    reads_SOA, mean_read_length, match_frequencies = analyse_dataset(config, ref_genome_file_path, fastq_file_paths)
    
    mkpath(output_dir)

    suffix = "k$(k)s$(step)"

    CSV.write("$output_dir/reads_result-$suffix.csv", reads_SOA)
    CSV.write("$output_dir/match_frequencies-$suffix.csv", match_frequencies)

    match_frequency_plot("$output_dir/match_frequencies-$suffix.csv", mean_read_length, step, "$output_dir/match_frequencies-$suffix.svg")

    #thresholds = [70, 100, 130]
    #activity_plot(
    #    generate_activity_vectors(
    #        reads_SOA,
    #        Int32(fasta_first_seqlen(ref_genome_file_path)),
    #        thresholds
    #    ),
    #    "$output_dir/activity-$suffix.png"
    #)
end

export snipe_reads