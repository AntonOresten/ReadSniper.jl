

function fetch_query_data(
    query_file::AbstractString,
    k::Int64,
)
    reader = FASTA.Reader(open(query_file, "r"))

    query_sequence = FASTA.sequence(LongDNA{4}, first(reader))
    query_length = length(query_sequence)
    query_kmer_vector = create_kmer_vector(query_sequence, k)

    close(reader)

    query_kmer_dict = kmer_index_dict(query_kmer_vector)

    return query_sequence, query_length, query_kmer_dict
end

export fetch_query_data


function fasta_first_seqlen(fasta_file::AbstractString)
    reader = FASTA.Reader(open(fasta_file, "r"))
    sequence = FASTA.sequence(LongDNA{4}, first(reader))
    len = length(sequence)
    close(reader)
    return len
end


function analyse_dataset(
    config::Config,
    query_file::AbstractString,
    datafiles_paths::Vector{<:AbstractString},
)
    query_sequence, query_length, query_kmer_dict = fetch_query_data(query_file, config.k)
    match_frequency_dict::Dict{Int32, Int32} = Dict()

    read_result_sets, mean_read_length = parallel_scanning(config, datafiles_paths, query_kmer_dict, match_frequency_dict)
    read_result_SOA = sort!(StructArray(vcat([collect(rr_set) for rr_set in read_result_sets]...)), by=result->result.kmer_matches, rev=true)
    # SOA: Structure-Of-Arrays

    sorted_match_bins = sort!(collect(match_frequency_dict), by=pair->pair[1])
    matches, frequency = getindex.(sorted_match_bins, 1), getindex.(sorted_match_bins, 2)

    return read_result_SOA, mean_read_length, (matches=matches, frequency=frequency)
end

export analyse_dataset


"""
Returns the reads whose indices are contained in a given set
"""
function filter_fastq(fastq_file, read_ids::Set{Int64})
    filename = split(fastq_file, ".")[1]
    output_file = filename*"-filtered.fasta"
    reader = FASTQ.Reader(open(fastq_file, "r"))
    writer = FASTA.Writer(open(output_file, "w"))
    for (read_id, record) in enumerate(reader)
        if read_id in read_ids
            write(writer, FASTA.Record(filename*".$read_id", FASTQ.sequence(LongDNA{4}, record)))
        end
    end
    close(writer)
    close(reader)
end

export filter_fastq


function count_fastq_records(fastq_file::AbstractString)
    k = 0
    record = FASTQ.Record()
    reader = FASTQ.Reader(open(fastq_file, "r"))
    while !eof(reader)
        read!(reader, record)
        k += 1
    end
    close(reader)
    return k
end

export count_fastq_records


#=
function quick_fastq_stats(fastq_file::AbstractString)
    reader = FASTQ.Reader(open(fastq_file_file, "r"))
    median_read_length = nothing
end
=#