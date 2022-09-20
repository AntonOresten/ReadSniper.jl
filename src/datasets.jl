
using CSV
using DataFrames

using FASTX
using BioSequences

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

function get_best_datasets(sra_metadata_file::String, n::Int64)
    datasets::DataFrame = DataFrames.sort(CSV.read(sra_metadata_file, DataFrame), [:pident], rev=true)
    run_ids = datasets.run_id
    best_SRR::Vector{String} = []
    for run_id in run_ids
        if length(best_SRR) >= n
            break
        end
        if occursin("SRR", run_id)
            push!(best_SRR, run_id)
        end
    end
    return best_SRR[1:n]
end

function download_datasets(datasets::Vector{String})
    nothing
end

#@show get_best_datasets("_tmp_submission.sra.metadata.csv", 20)
