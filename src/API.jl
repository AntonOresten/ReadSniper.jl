

function snipe_reads(
    ref_genome_file_path::String,
    fastq_file_dir::String,
    fastq_files::Vector{String};
    output_dir::String="",
    k::Int64=7,
    step::Int64=1,
)
    read_scores_SOA, score_counts = score_dataset(ref_genome_file_path, fastq_file_dir, fastq_files, k=k, step=step)

    create_dir(output_dir)

    suffix = "k$(k)s$(step)"

    CSV.write("$output_dir/read_scores-$suffix.csv", read_scores_SOA)
    CSV.write("$output_dir/score_counts-$suffix.csv", score_counts)

    score_frequency_plot(score_counts.scores, score_counts.counts, "$output_dir/score_frequencies-$suffix.png")

    thresholds = [70, 100, 130]
    activity_plot(
        generate_activity_vectors(
            read_scores_SOA,
            Int32(fasta_first_seqlen(ref_genome_file_path)),
            thresholds
        ),
        "$output_dir/activity-$suffix.png"
    )
end

export snipe_reads