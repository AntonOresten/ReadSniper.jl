

function fetch_query_data(
    query_file::String,
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


function fasta_first_seqlen(fasta_file::String)
    reader = FASTA.Reader(open(fasta_file, "r"))
    sequence = FASTA.sequence(LongDNA{4}, first(reader))
    len = length(sequence)
    close(reader)
    return len
end


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


function count_fastq_records(fastq_file::String)
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

export fqRecordCount


function split_read_sequence(seq::LongDNA{4})
    half_len = length(seq) ÷ 2
    return seq[1:half_len], seq[half_len+1:end]
end

# TODO: Function for getting the best datasets out of a CSV file with pident and coverage values?

#=
function quick_fastq_stats(fastq_file::String)
    reader = FASTQ.Reader(open(fastq_file_file, "r"))
    median_read_length = nothing
end
=#