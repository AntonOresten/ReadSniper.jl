
module ReadSniper


include("include_all.jl")

import Dates
using ProgressBars


function retrieve_reads(
    reference_file::String,
    SRR_ID::String,
    split_files_path::String,
    score_threshold::Float64=0.10,
    k::Int64=10,
    step::Int64=1)

    output_dir = "$SRR_ID-$(Dates.now())"

    reference_reader = FASTA.Reader(open(reference_file, "r"))
    reference_sequence = FASTA.sequence(LongDNASeq, first(reference_reader))
    reference_length = length(reference_sequence)

    reference_kmer_vector = create_kmer_vector(reference_sequence, k)
    reference_kmer_dict = single_match_kmer_dict(reference_kmer_vector)

    close(reference_reader)

    file1 = SRR_ID*"_1"
    file2 = SRR_ID*"_2"

    reader1 = FASTQ.Reader("split_files_path\\$file1.fastq", "r")
    reader2 = FASTQ.Reader("split_files_path\\$file2.fastq", "r")

    writer = FASTA.Writer("$output_dir\\$SRR_ID-filtered.fastq", "w")

    for (read_id, (record1::FASTQ.Record, record2::FASTQ.Record)) in ProgressBar(enumerate(zip(reader1, reader2)))
        read1_length = FASTQ.seqlen(record1)
        read2_length = FASTQ.seqlen(record2)

        read1::LongDNASeq = FASTQ.sequence(LongDNASeq, record1)
        read1_match_indices::Vector{Int64} = get_single_match_indices(read1, reference_kmer_dict, k, step)
        read1_score::Float64 = max_window(read1_match_indices, read1_length) / read1_length

        read2::LongDNASeq = FASTQ.sequence(LongDNASeq, record2)
        read2_match_indices::Vector{Int64} = get_single_match_indices(read2, reference_kmer_dict, k, step)
        read2_score::Float64 = max_window(read2_match_indices, read2_length) / read2_length

        if read1_score > score_threshold || read2_score > score_threshold
            write(writer, FASTA.Record("$(file1)_$read_id*", read1))
            write(writer, FASTA.Record("$(file2)_$read_id*", read2))
        end
    end
    number_of_reads = read_id

    close(writer)

    close(reader1)
    close(reader2)
    
end

retrieve_reads(
    "reference.fasta",
    "SRR10873757",
    "SRR10873757",
    0.1,
)

export

# utils.jl

size_score,
maximum_increasing_subsequence,
max_window,    

# kmers.jl

create_kmer_vector,
single_match_kmer_dict,
get_single_match_indices


end # module