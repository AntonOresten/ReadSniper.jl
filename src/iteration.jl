

function parallel_scanning(
    config::Config,
    reference::Reference,
    fasta_metadata_list::Vector{FASTAMetadata},
    match_frequency_dict::Dict{Int32, Int32},
    show_progress::Bool=true,
)
    #false_vector = zeros(Bool, reference.length-config.k+1)
    threshold = config.threshold

    file_count = length(fasta_metadata_list)
    read_count = sum([fasta_metadata.read_count for fasta_metadata in fasta_metadata_list])
    threads_utilized = min(config.nthreads, file_count)

    reads_scanned = zeros(Int, file_count)

    read_result_sets = [Set{ReadResult}() for _ in 1:file_count]

    show_progress ? progress_bar = ProgressBar(config, read_count, file_count, threads_utilized) : nothing
    Threads.@threads for (i, fasta_metadata) in collect(enumerate(fasta_metadata_list))
        read_number::Int = 0
        reader = FASTAReader(open(fasta_metadata.path, "r"); copy=false)
        while !eof(reader)
            record = first(reader)
            reads_scanned[i] = (read_number += 1)

            read_result = scan_read(config, read_number, record, reference.kmer_dict)

            increment_dict_value!(read_result.kmer_matches, match_frequency_dict)

            if read_result.kmer_matches > threshold
                push!(read_result_sets[i], read_result)
            end

            if show_progress && (i == 1 || threads_utilized == 1) && iszero(read_number % 1729)
                update(progress_bar, sum(reads_scanned), sum(length(reads) for reads in read_result_sets))
            end
        end
        close(reader)
    end

    show_progress ? finish(progress_bar, sum(length(reads) for reads in read_result_sets)) : nothing

    return read_result_sets
end


# scattered scanning


function analyse_dataset(
    config::Config,
    reference::Reference,
    fasta_metadata_list::Vector{FASTAMetadata},
)
    match_frequency_dict = Dict{Int32, Int32}()

    read_result_sets = parallel_scanning(config, reference, fasta_metadata_list, match_frequency_dict)
    read_result_SOA = sort!(StructArray(vcat([collect(rr_set) for rr_set in read_result_sets]...)), by=result->result.kmer_matches, rev=true)

    sorted_match_bins = sort!(collect(match_frequency_dict), by=pair->pair[1])
    matches, frequency = getindex.(sorted_match_bins, 1), getindex.(sorted_match_bins, 2)

    return read_result_SOA, (matches=matches, frequency=frequency)
end

export analyse_dataset
