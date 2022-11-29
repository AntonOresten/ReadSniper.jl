
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
            write(writer, FASTA.Record(filename*".$read_id", FASTQ.sequence(LongDNASeq, record)))
        end
    end
    close(writer)
    close(reader)
end

function fqRecordCount(fastq_file::String)
    reader = FASTQ.Reader(open(fastq_file, "r"))
    k = 0
    for record in reader
        k += 1
    end
    close(reader)
    return k
end

function faRecordCount(fasta_file::String)
    nothing
end

# TODO: Function for getting the best datasets out of a CSV file with pident and coverage values?