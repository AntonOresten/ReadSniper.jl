

function parallel_scanning(
    config::Config,
    reference::Reference,
    datafile_list::Vector{DatafileMetadata},
    score_frequency_dict::Dict{Int32, Int32},
    show_progress::Bool = true,
)
    file_count = length(datafile_list)
    read_count = sum([datafile.read_count for datafile in datafile_list])
    reads_scanned = zeros(Int, file_count)
    threads_utilized = min(config.nthreads, file_count)

    read_result_sets = [Set{ReadResult}() for _ in 1:file_count]

    if show_progress progress_bar = ProgressBar(config, read_count, file_count, threads_utilized) end
    Threads.@threads for (i, datafile) in collect(enumerate(datafile_list))
        read_number = 0
        reader = FASTAReader(open(datafile.path, "r"); copy=false)
        while !eof(reader)
            record = first(reader)
            reads_scanned[i] = (read_number += 1)

            ref_range_start, ref_range_end, score = scan_read(config, record, reference.kmer_dict)
            increment_dict_value!(score, score_frequency_dict)
            if score >= config.threshold
                push!(read_result_sets[i], ReadResult(read_number, ref_range_start, ref_range_end, score))
            end

            if show_progress && (i == 1 || threads_utilized == 1) && iszero(read_number % 1729)
                update(progress_bar, sum(reads_scanned), sum(length(reads) for reads in read_result_sets))
            end
        end
        close(reader)
    end

    if show_progress finish(progress_bar, sum(length(reads) for reads in read_result_sets)) end

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