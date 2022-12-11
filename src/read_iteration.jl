
# TODO: add struct for Reader with given start and end

function progress_bar_values(
    reads_scanned::Int64, reads_total::Int64, reads_sniped::Int64,
    file_number::Int64, number_of_files::Int64,
    k::Int64, step::Int64, match_method::Type, number_of_threads::Int64,
)
    return [
        (Symbol("Reads scanned"), "$reads_scanned/$reads_total"), (Symbol("Reads sniped"), reads_sniped),
        (Symbol("Files scanned"), "$file_number/$number_of_files"),
        (Symbol("k-mer length"), k), (Symbol("k-mer step"), step),
        (Symbol("Match method"), (match_method == Int64 ? "Single" : "Multi")), (Symbol("Threads"), number_of_threads),
    ]
end


function scan_dataset(
    dataset::Vector{String},
    query_kmer_dict::Dict{LongSubSeq{DNAAlphabet{4}}, T},
    k::Int64,
    score_distribution_dict::Dict{Int64, Int64};
    step::Int64,
    threshold::Int64,
    number_of_threads::Int64=1,
) where T <: Union{Int64, Vector{Int64}}

    k_1 = k - 1
    threshold = ceil(Int64, threshold / step)

    record_count = fqRecordCount(dataset[1])
    datafile_count = length(dataset)
    read_count = record_count * datafile_count
    
    # Preventing threads from adding to one integer/set at the same time
    reads_scanned = zeros(Int64, number_of_threads)
    reads_sniped = zeros(Int64, number_of_threads)
    reads = [Set{ReadMetadata}() for _ in 1:number_of_threads]

    window = record_count ÷ number_of_threads
    remainder = record_count % number_of_threads

    println()
    progress_bar = Progress(read_count, desc="Sniping reads... ", color=:white)
    for (datafile_num, datafile) in enumerate(dataset)

        # DatafileReaders each have a designated window of each dataset
        datafile_readers::Vector{DatafileReader} = [
            DatafileReader(FASTQ.Reader(open(datafile, "r")), reader_id, window, remainder)
            for reader_id in 1:number_of_threads]
        
        #for df in datafile_readers
        #    println(df)
        #end

        @sync begin
            Threads.@threads for reader_id in 1:number_of_threads
                df_reader = datafile_readers[reader_id]
                window = df_reader.window
                reader = df_reader.reader

                for (read_number, record) in enumerate(reader)
                    if read_number > window break end

                    max_matches_in_range, range_start, range_end, query_indices = analyze_record(record, query_kmer_dict, k_1, step, threshold)

                    reads_scanned[reader_id] += 1
                    increment_dict_value!(max_matches_in_range, score_distribution_dict)

                    if max_matches_in_range > threshold
                        read = ReadMetadata(df_reader.first_index+read_number-1, datafile_num, step*max_matches_in_range, query_indices[range_start], query_indices[range_end])
                        push!(reads[reader_id], read)
                        reads_sniped[reader_id] += 1
                    end

                    if (reader_id == 1) && iszero(read_number % 369)
                        sum_reads_scanned = sum(reads_scanned)
                        update!(progress_bar, sum_reads_scanned, showvalues = progress_bar_values(
                            sum_reads_scanned, read_count, sum(reads_sniped),
                            datafile_num, datafile_count,
                            k, step, T, number_of_threads))
                    end
                end
            end
        end
    end
    
    update!(progress_bar, read_count, :green, showvalues = progress_bar_values(
        sum(reads_scanned), read_count, sum(reads_sniped),
        datafile_count, datafile_count,
        k, step, T, number_of_threads))
    finish!(progress_bar)
    return reads
end


function retrieve_reads(
    query_file::String,
    dataset::Vector{String};
    output_dir::String="output",
    k::Int64=10,
    step::Int64=1,
    single_match::Bool=false,
)
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    number_of_threads = Threads.nthreads()
    #addprocs(number_of_threads)

    query_sequence, query_length, query_kmer_dict = fetch_query_data(query_file, k, single_match)
    score_distribution_dict::Dict{Int64, Int64} = Dict()

    read_set_vector = scan_dataset(
        dataset, query_kmer_dict, k, score_distribution_dict,
        step=step, threshold=9, number_of_threads=number_of_threads
    )
    read_struct_array = StructArray(vcat([collect(read_set) for read_set in read_set_vector]...))

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