
using BioSequences
using FASTX


module reads

include("utils.jl")
include("kmers.jl")

import Dates
using ProgressBars


function retrieve_reads_sm(
        reference_file::String,
        SRR_ID::String,
        split_files_path::String,
        score_threshold::Float64=0.10,
        k::Int64=10,
        step::Int64=1)
    """
    Uses a single-match dictionary when retrieving reads
    """

    output_dir = "$SRR_ID-$(Dates.today())"

    reference_reader = FASTA.Reader(open(reference_file, "r"))
    reference_sequence = FASTA.sequence(LongDNASeq, first(reference_reader))
    reference_length = length(reference_sequence)

    reference_kmer_vector = create_kmer_vector(reference_sequence, k)
    reference_kmer_dict = single_match_kmer_dict(reference_kmer_vector)

    close(reference_reader)

    file1 = SRR_ID*"_1"
    file2 = SRR_ID*"_2"

    reader1 = FASTQ.Reader(open("$split_files_path\\$file1.fastq", "r"))
    reader2 = FASTQ.Reader(open("$split_files_path\\$file2.fastq", "r"))

    if !isdir(output_dir)
        mkdir(output_dir)
    end
    writer = FASTA.Writer(open("$output_dir\\$SRR_ID-filtered.fasta", "w"))

    for (read_number, (record1::FASTQ.Record, record2::FASTQ.Record)) in ProgressBar(enumerate(zip(reader1, reader2)))
        read1_length = FASTQ.seqlen(record1)
        read2_length = FASTQ.seqlen(record2)

        read1::LongDNASeq = FASTQ.sequence(LongDNASeq, record1)
        read1_match_indices::Vector{Int64} = get_single_match_indices(read1, reference_kmer_dict, step)
        read1_score::Float64 = max_window(read1_match_indices, read1_length) / read1_length

        read2::LongDNASeq = FASTQ.sequence(LongDNASeq, record2)
        read2_match_indices::Vector{Int64} = get_single_match_indices(read2, reference_kmer_dict, step)
        read2_score::Float64 = max_window(read2_match_indices, read2_length) / read2_length

        if read1_score > score_threshold || read2_score > score_threshold
            write(writer, FASTA.Record("$(SRR_ID)_1_$read_number", read1))
            write(writer, FASTA.Record("$(SRR_ID)_2_$read_number", read2))
        end
    end
    close(writer)
    close(reader1)
    close(reader2)
end

#=retrieve_reads_sm(
    "queries/picorna_nuevo_rc.fasta",
    "SRR10873757",
    "SRR10873757",
    0.1)=#

function retrieve_reads_mm(
        reference_file::String,
        SRR_ID::String,
        split_files_path::String,
        score_threshold::Float64=0.10,
        k::Int64=10,
        step::Int64=1)
    """
    Uses a multi-match dictionary when retrieving reads
    """
    output_dir = "$SRR_ID-$(Dates.today())"

    reference_reader = FASTA.Reader(open(reference_file, "r"))
    reference_sequence = FASTA.sequence(LongDNASeq, first(reference_reader))
    reference_length = length(reference_sequence)

    reference_kmer_vector = create_kmer_vector(reference_sequence, k)
    reference_kmer_dict = multi_match_kmer_dict(reference_kmer_vector)

    close(reference_reader)

    file1 = SRR_ID*"_1"
    file2 = SRR_ID*"_2"

    reader1 = FASTQ.Reader(open("$split_files_path\\$file1.fastq", "r"))
    reader2 = FASTQ.Reader(open("$split_files_path\\$file2.fastq", "r"))

    if !isdir(output_dir)
        mkdir(output_dir)
    end
    writer = FASTA.Writer(open("$output_dir\\$SRR_ID-filtered.fastq", "w"))

    for (read_number, (record1::FASTQ.Record, record2::FASTQ.Record)) in ProgressBar(enumerate(zip(reader1, reader2)))
        read1_length = FASTQ.seqlen(record1)
        read2_length = FASTQ.seqlen(record2)

        read1::LongDNASeq = FASTQ.sequence(LongDNASeq, record1)
        read1_match_indices::Vector{Vector{Int64}} = get_multi_match_indices(read1, reference_kmer_dict, k, step)
        read1_score::Float64 = max_window(read1_match_indices, read1_length) / read1_length

        read2::LongDNASeq = FASTQ.sequence(LongDNASeq, record2)
        read2_match_indices::Vector{Vector{Int64}} = get_multi_match_indices(read2, reference_kmer_dict, k, step)
        read2_score::Float64 = max_window(read2_match_indices, read2_length) / read2_length

        if read1_score > score_threshold || read2_score > score_threshold
            write(writer, FASTA.Record("$(SRR_ID)_1_$read_number", read1))
            write(writer, FASTA.Record("$(SRR_ID)_2_$read_number", read2))
        end
    end
    close(writer)
    close(reader1)
    close(reader2)
end

export

retrieve_reads_sm,
retrieve_reads_mm

end