

function parallel_scanning(
    config::Config,
    reference::Reference,
    datafile_list::Vector{DatafileMetadata},
    score_frequency_dict_union::Dict{Int32, Int32},
    show_progress::Bool = true,
)
    file_count = length(datafile_list)
    read_count = sum([datafile.read_count for datafile in datafile_list])
    reads_scanned = zeros(Int, config.nthreads)

    read_result_sets = [Set{ReadResult}() for _ in 1:config.nthreads]
    score_frequency_dicts = [Dict{Int32, Int32}() for _ in 1:config.nthreads]

    if show_progress progress_bar = ProgressBar(config, file_count, read_count, config.nthreads) end
    for (file_number, datafile) in collect(enumerate(datafile_list))
        if show_progress update(progress_bar, file_number-1, sum(reads_scanned), sum(length(rs) for rs in read_result_sets)) end

        idx = faidx(open(datafile.path, "r"))
        readers = [FASTAReader(open(datafile.path, "r"), index=idx, copy=false) for _ in 1:config.nthreads]
        number_of_readers = length(readers)
        reader_window_size = datafile.read_count ÷ number_of_readers

        Threads.@threads for (i, reader) in collect(enumerate(readers))
            window_start = (i - 1) * reader_window_size + 1
            window_end = i * reader_window_size + (i == number_of_readers ? datafile.read_count % number_of_readers : 0)
            seekrecord(reader, window_start)
            for read_number in window_start:window_end
                record = first(reader)
                reads_scanned[i] += 1

                ref_range_start, ref_range_end, score = scan_read(config, record, reference.kmer_dict)
                increment_dict_value!(score, 1, score_frequency_dicts[i])
                if score >= config.threshold
                    push!(read_result_sets[i], ReadResult(read_number, ref_range_start, ref_range_end, score))
                end

                if show_progress && i == 1 && iszero(read_number % 729)
                    update(progress_bar, file_number-1, sum(reads_scanned), sum(length(rs) for rs in read_result_sets))
                end
            end
        end
        for reader in readers close(reader) end
        idx = nothing
        GC.gc()
    end

    for sfd in score_frequency_dicts
        for (k, v) in sfd
            increment_dict_value!(k, v, score_frequency_dict_union)
        end
    end

    if show_progress finish(progress_bar, sum(length(rs) for rs in read_result_sets)) end

    return read_result_sets
end


function analyse_dataset(
    config::Config,
    reference::Reference,
    datafile_list::Vector{DatafileMetadata},
    show_progress::Bool = true,
)
    score_frequency_dict = Dict{Int32, Int32}()

    read_result_sets = parallel_scanning(config, reference, datafile_list, score_frequency_dict, show_progress)
    read_result_SOA = sort!(StructArray(vcat([collect(rr_set) for rr_set in read_result_sets]...)), by=result->result.score, rev=true)

    sorted_match_bins = sort!(collect(score_frequency_dict), by=pair->pair[1])
    matches, frequency = getindex.(sorted_match_bins, 1), getindex.(sorted_match_bins, 2)

    return read_result_SOA, (matches=matches, frequency=frequency)
end

export analyse_dataset


# write all read scores to a file maybe - 16 bits per score, and then enumerate when reading to get read num