
# TODO: add struct for Reader with given start and end

function progress_bar_values(
    reads_sniped::Int64, reads_scanned::Int64, reads_total::Int64,
    k::Int64, step::Int64, match_method::Type, threads_used::Int64,
)
    return [
        (Symbol("Reads sniped"), reads_sniped), (Symbol("Reads scanned"), reads_scanned),
        (Symbol("Number of reads"), reads_total), (Symbol("k-mer length"), k), (Symbol("k-mer step"), step),
        (Symbol("Match method"), (match_method == Int64 ? "Single" : "Multi")), (Symbol("Threads"), threads_used),
    ]
end

function scan_dataset(
    readers::Vector, 
    record_count::Int64,
    query_kmer_dict::Dict{LongSubSeq{DNAAlphabet{4}}, T},
    k::Int64,
    score_distribution_dict::Dict{Int64, Int64};
    step::Int64=1,
    threshold::Int64=10,
) where T <: Union{Int64, Vector{Int64}}

    k_ = k - 1
    threshold = ceil(Int64, threshold / step)

    reader_count = length(readers)
    read_count = record_count * reader_count
    reads::Vector{Set{ReadMetadata}} = [Set() for _ in readers]
    
    # Vector{Int} instead of Int to prevent two readers from adding to one at the same time
    reads_sniped = zeros(Int64, reader_count)

    println()
    threads_used = Threads.nthreads()
    progress_bar = Progress(read_count, desc="Sniping reads... ", color=:white)
    Threads.@threads for (i, reader) in collect(enumerate(readers))
        for read_number in 1:record_count
            record = iterate(reader)[1]

            max_matches_in_range, range_start, range_end, query_indices = analyze_record(record, query_kmer_dict, k_, step, threshold)
            
            if max_matches_in_range > threshold
                read = ReadMetadata(read_number, i, max_matches_in_range, query_indices[range_start], query_indices[range_end])
                push!(reads[i], read)
                reads_sniped[i] += 1
            end

            increment_dict_value!(max_matches_in_range, score_distribution_dict)
            
            if (i == 1) && (read_number % 369 == 0)
                update!(progress_bar, 2*read_number, showvalues = progress_bar_values(
                    sum(reads_sniped), 2*read_number, read_count,
                    k, step, T, threads_used))
            end
            #push!(read_ids[i], read_num)
            #stack_index_range!(coverage, sorted_indices[start_index]:sorted_indices[end_index])
        end
    end

    update!(progress_bar, read_count, :green, showvalues = progress_bar_values(
        sum(reads_sniped), read_count, read_count,
        k, step, T, threads_used))
    finish!(progress_bar)
    println()

    return reads
end


function retrieve_reads(
    query_file::String,
    dataset::Tuple{String, String};
    output_dir::String="output",
    k::Int64=10,
    step::Int64=1,
    single_match::Bool=false,
)
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    query_sequence, query_length, query_kmer_dict = fetch_query_data(query_file, k, single_match)
    score_distribution_dict::Dict{Int64, Int64} = Dict()

    # Two sets of reads due to paired-end sequencing
    dataset_1, dataset_2 = dataset

    reader1 = FASTQ.Reader(open("$dataset_1", "r"))
    reader2 = FASTQ.Reader(open("$dataset_2", "r"))

    record_count = fqRecordCount(dataset_1)

    reads_1, reads_2 = scan_dataset([reader1, reader2], record_count, query_kmer_dict, k, score_distribution_dict, step=step, threshold=10)
    read_struct_array = StructArray([reads_1..., reads_2...])

    println("Saving results...")
    CSV.write("$output_dir/reads.csv", read_struct_array)

    println("Saving read match distribution...")
    sorted_score_bins = sort(collect(score_distribution_dict), by=x->x[1])
    
    scores = collect([pair[1] for pair in sorted_score_bins])
    counts = collect([pair[2] for pair in sorted_score_bins])

    CSV.write("$output_dir/read_match_distribution.csv", (scores=scores, counts=counts))

    #println("Generating activity vectors...")
    #activity_vectors = generate_activity_vectors(read_struct_array, query_length, [10,20,30])
    #println("Creating histogram...")
    #plot_activity_histogram(activity_vectors)

    #bar_plot(scores, counts, "read_match_distribution.pdf")

    #CSV.write("$output_dir/read_ids.csv", DataFrame(read_ids=sort(collect(union(read_ids_1, read_ids_2)))))
    #CSV.write("$output_dir/coverage.csv", DataFrame(i=1:query_length, frequency=coverage))
end

export retrieve_reads