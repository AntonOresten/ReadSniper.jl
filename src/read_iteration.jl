

function parallel_scanning(
    config::Config,
    datafiles_paths::Vector{String},
    kmer_index_dict::Dict{LongSubSeq{DNAAlphabet{4}}, Vector{Int64}},
    score_distribution_dict::Dict{Int32, Int32},
)

    file_count = length(datafiles_paths)
    read_count = file_count * count_fastq_records(datafiles_paths[1])
    threads_utilized = min(config.nthreads, file_count)

    reads_scanned = zeros(Int32, file_count)
    read_result_sets = [Set{ReadResult}() for _ in 1:file_count]

    progress_bar = ProgressBar(config, read_count, file_count, threads_utilized)
    Threads.@threads for (i, datafile) in collect(enumerate(datafiles_paths))
        read_number::Int32 = 0
        reader = FASTQ.Reader(open(datafile, "r"); copy=false)
        while !eof(reader)
            reads_scanned[i] = (read_number += 1)
            record = first(reader)

            read_result = scan_read(config, read_number, record, kmer_index_dict)
            increment_dict_value!(read_result.kmer_match_count, score_distribution_dict)

            if read_result.kmer_match_count > config.threshold
                push!(read_result_sets[i], read_result)
            end

            if (i == 1) && iszero(read_number % 1729)
                update(progress_bar, sum(reads_scanned))
            end
        end
        close(reader)
    end
    finish(progress_bar)

    return read_result_sets
end