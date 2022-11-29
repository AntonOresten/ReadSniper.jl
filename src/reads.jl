

function read_query(
        query_file::String,
        k::Int,
        single_match::Bool=false)

    query_reader = FASTA.Reader(open(query_file, "r"))

    query_sequence = FASTA.sequence(LongDNASeq, first(query_reader))
    query_length = length(query_sequence)
    query_kmer_vector = create_kmer_vector(query_sequence, k)

    close(query_reader)

    if single_match
        query_kmer_dict = single_match_kmer_dict(query_kmer_vector)
    else
        query_kmer_dict = multi_match_kmer_dict(query_kmer_vector)
    end

    return query_sequence, query_length, query_kmer_dict
end


function split_read_sequence(seq::LongDNASeq)
    half_len = length(seq) ÷ 2
    return seq[1:half_len], seq[half_len+1:end]
end


"""
Uses a single-match dictionary when retrieving reads.

Works best with small genomes (or large k-values) due to k-mers only being assigned one index.
"""
function retrieve_reads_sm(
    query_file::String,
    dataset::Union{String, Tuple{String, String}},
    score_threshold::Float64=0.10,
    k::Int64=10,
    step::Int64=1,
)

    output_dir = "$dataset-$(Dates.today())"

    file1 = dataset*"_1.fastq"
    file2 = dataset*"_2.fastq"

    reader1 = FASTQ.Reader(open("$path_to_dataset\\$file1", "r"))
    reader2 = FASTQ.Reader(open("$path_to_dataset\\$file2", "r"))

    if !isdir(output_dir)
        mkdir(output_dir)
    end
    writer = FASTA.Writer(open("$output_dir\\$dataset-filtered.fasta", "w"))

    for (read_number, (record1::FASTQ.Record, record2::FASTQ.Record)) in ProgressBar(enumerate(zip(reader1, reader2)))
        read1_length = FASTQ.seqlen(record1)
        read2_length = FASTQ.seqlen(record2)

        read1::LongDNASeq = FASTQ.sequence(LongDNASeq, record1)
        read1_match_indices::Vector{Int64} = single_match_indices(read1, query_kmer_dict, step)
        read1_score::Float64 = max_window(read1_match_indices, read1_length) / read1_length

        read2::LongDNASeq = FASTQ.sequence(LongDNASeq, record2)
        read2_match_indices::Vector{Int64} = single_match_indices(read2, query_kmer_dict, step)
        read2_score::Float64 = max_window(read2_match_indices, read2_length) / read2_length

        if read1_score > score_threshold || read2_score > score_threshold
            write(writer, FASTA.Record("$(dataset_name)_1_$read_number", read1))
            write(writer, FASTA.Record("$(dataset_name)_2_$read_number", read2))
        end
    end
    close(writer)
    close(reader1)
    close(reader2)
end


"Single-match scanning"
function scan_reads(
    reader1::FASTQ.Reader,
    reader2::FASTQ.Reader,
    query_kmer_dict::Dict{DNAMer{k}, Int64},
    coverage::Vector{Int64},
    score_distribution_dict::Dict{Int64, Int64},
    step::Int64=1,
    threshold::Int64=25,
) where {k}

    read_ids::Vector{Set{Int64}} = [Set(), Set()]

    # i will be either 1 or 2
    Threads.@threads for (i, reader) in collect(enumerate((reader1, reader2)))
        for (read_num, record::FASTQ.Record) in ProgressBar(enumerate(reader))
            
            read_length = FASTQ.seqlen(record)
            read::LongDNASeq = FASTQ.sequence(LongDNASeq, record)
            read_match_indices::Vector{Int64} = single_match_indices(read, query_kmer_dict, step)
            read_score::Int64, start_index, end_index = max_window(read_match_indices, read_length, lis)
            increment_value(score_distribution_dict, read_score)
            
            if read_score ≥ threshold
                push!(read_ids[i], read_num)
                stack_index_range!(coverage, read_match_indices[start_index]:read_match_indices[end_index])
            end
        end
    end

    return read_ids
end


struct Read
    id::Int64
    matches::Int64
    match_range_start::Int64
    match_range_end::Int64
end


"Multi-match scanning"
function scan_reads(
    reader1::FASTQ.Reader,
    reader2::FASTQ.Reader,
    query_kmer_dict::Dict{DNAMer{k}, Vector{Int64}},
    coverage::Vector{Int64},
    score_distribution_dict::Dict{Int64, Int64},
    step::Int64=1,
    threshold::Int64=25,
) where {k}

    read_ids::Vector{Set{Int64}} = [Set(), Set()]

    # i will be either 1 or 2
    Threads.@threads for (i, reader) in collect(enumerate((reader1, reader2)))
        for (read_num, record::FASTQ.Record) in ProgressBar(enumerate(reader))
            
            read_length = FASTQ.seqlen(record)
            read::LongDNASeq = FASTQ.sequence(LongDNASeq, record)
            read_match_index_vectors::Vector{Vector{Int64}} = multi_match_indices(read, query_kmer_dict, step)
            concat_indices::Vector{Int64} = vcat(read_match_index_vectors...)
            sorted_indices::Vector{Int64} = sort(concat_indices)
            matches_in_window::Int64, start_index::Int64, end_index::Int64 = max_window(sorted_indices, read_length, length)
            increment_value(score_distribution_dict, matches_in_window)
            
            push!(read_ids[i], read_num)
            stack_index_range!(coverage, sorted_indices[start_index]:sorted_indices[end_index])
        end
    end

    return read_ids
end




function retrieve_reads(
    query_file::String,
    dataset::Tuple{String, String};
    output_dir::String="output",
    k::Int64=10,
    threshold::Int64=25,
    step::Int64=1,
    single_match::Bool=false,
)
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    query_sequence, query_length, query_kmer_dict = read_query(query_file, k, single_match)
    coverage = zeros(Int64, query_length)
    score_distribution_dict::Dict{Int64, Int64} = Dict()

    # Two sets of reads because of paired-end sequencing
    reads_1, reads_2 = dataset

    read_count = fqRecordCount(reads_1)

    reader1 = FASTQ.Reader(open("$reads_1", "r"))
    reader2 = FASTQ.Reader(open("$reads_2", "r"))

    output_filename = split(reads_1, '_')[1]*".fastq"

    read_ids_1, read_ids_2 = scan_reads(reader1, reader2, query_kmer_dict, coverage, score_distribution_dict, step, threshold)

    CSV.write("$output_dir/read_ids.csv", DataFrame(read_ids=sort(collect(union(read_ids_1, read_ids_2)))))
    CSV.write("$output_dir/coverage.csv", DataFrame(i=1:query_length, frequency=coverage))
end

export retrieve_reads