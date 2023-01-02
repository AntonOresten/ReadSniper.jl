

function parallel_scanning(
    config::Config,
    datafiles_paths::Vector{<:AbstractString},
    kmer_index_dict::Dict{LongSubSeq{DNAAlphabet{4}}, Vector{Int64}},
    match_frequency_dict::Dict{Int32, Int32},
)
    file_count = length(datafiles_paths)
    read_count = file_count * count_fastq_records(datafiles_paths[1])
    threads_utilized = min(config.nthreads, file_count)

    reads_scanned = zeros(Int32, file_count)
    bases_parsed = zeros(Int32, file_count)

    read_result_sets = [Set{ReadResult}() for _ in 1:file_count]

    progress_bar = ProgressBar(config, read_count, file_count, threads_utilized)
    Threads.@threads for (i, datafile) in collect(enumerate(datafiles_paths))
        read_number::Int = 0
        reader = FASTQ.Reader(open(datafile, "r"))
        while !eof(reader)
            record = first(reader)

            reads_scanned[i] = (read_number += 1)
            bases_parsed[i] += FASTQ.seqsize(record)

            read_result = scan_read(config, read_number, record, kmer_index_dict)
            increment_dict_value!(read_result.kmer_matches, match_frequency_dict)

            if read_result.kmer_matches > config.threshold
                push!(read_result_sets[i], read_result)
            end

            if (i == 1) && iszero(read_number % 1729)
                update(progress_bar, sum(reads_scanned), sum(bases_parsed), sum(bases_parsed)÷sum(reads_scanned))
            end
        end
        close(reader)
    end
    mean_read_length = sum(bases_parsed) ÷ sum(reads_scanned)

    finish(progress_bar, sum(bases_parsed), mean_read_length)

    return read_result_sets, mean_read_length
end