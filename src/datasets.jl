

function fetch_query_data(
    query_file::String,
    k::Int64,
    single_match::Bool=false,
)

    reader = FASTA.Reader(open(query_file, "r"))

    query_sequence = FASTA.sequence(LongDNA{4}, first(reader))
    query_length = length(query_sequence)
    query_kmer_vector = create_kmer_vector(query_sequence, k)

    close(reader)

    if single_match
        query_kmer_dict = single_match_kmer_dict(query_kmer_vector)
    else
        query_kmer_dict = multi_match_kmer_dict(query_kmer_vector)
    end

    return query_sequence, query_length, query_kmer_dict
end

export fetch_query_data


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

function fqRecordCount(fastq_file::String)
    reader = FASTQ.Reader(open(fastq_file, "r"))
    k = 0
    for _ in reader
        k += 1
    end
    close(reader)
    return k
end

# TODO: Function for getting the best datasets out of a CSV file with pident and coverage values?