

function read_query(
        query_file::String,
        k::Int,
        single_match::Bool=false)

    query_reader = FASTA.Reader(open(query_file, "r"))

    query_sequence = FASTA.sequence(LongDNASeq, first(query_reader))
    query_length = length(query_sequence)

    query_kmer_vector = create_kmer_vector(query_sequence, k)

    if single_match
        query_kmer_dict = single_match_kmer_dict(query_kmer_vector)
    else
        query_kmer_dict = multi_match_kmer_dict(query_kmer_vector)
    end
    
    close(query_reader)

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


"""
Uses a multi-match dictionary when retrieving reads
"""
function retrieve_reads_mm(
    query_file::String,
    dataset::String,
    score_threshold::Float64=0.10,
    k::Int64=10,
    step::Int64=1,
)

    dataset_name = splitext(dataset)
    output_dir = dataset_name*"-$(Dates.today())"

    query_reader = FASTA.Reader(open(query_file, "r"))
    query_sequence = FASTA.sequence(LongDNASeq, first(query_reader))
    query_length = length(query_sequence)

    query_kmer_vector = create_kmer_vector(query_sequence, k)
    query_kmer_dict = multi_match_kmer_dict(query_kmer_vector)

    close(query_reader)

    file1 = dataset_name*"_1"
    file2 = dataset_name*"_2"

    reader1 = FASTQ.Reader(open("$split_files_path\\$file1.fastq", "r"))
    reader2 = FASTQ.Reader(open("$split_files_path\\$file2.fastq", "r"))

    if !isdir(output_dir)
        mkdir(output_dir)
    end
    writer = FASTA.Writer(open("$output_dir\\$dataset_name-filtered.fasta", "w"))

    for (read_number, (record1::FASTQ.Record, record2::FASTQ.Record)) in ProgressBar(enumerate(zip(reader1, reader2)))
        read1_length = FASTQ.seqlen(record1)
        read2_length = FASTQ.seqlen(record2)

        read1::LongDNASeq = FASTQ.sequence(LongDNASeq, record1)
        read1_match_indices::Vector{Int64} = multi_match_indices(read1, query_kmer_dict, step)
        read1_score::Float64 = max_window(read1_match_indices, read1_length) / read1_length

        read2::LongDNASeq = FASTQ.sequence(LongDNASeq, record2)
        read2_match_indices::Vector{Int64} = multi_match_indices(read2, query_kmer_dict, step)
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
    step::Int64=1,
    threshold::Int64=25,
) where {k}

    dataset1_read_numbers::Set{Int64} = Set([])
    dataset2_read_numbers::Set{Int64} = Set([])

    for (read_num, (record1::FASTQ.Record, record2::FASTQ.Record)) in ProgressBar(enumerate(zip(reader1, reader2)))
        
        read1_length = FASTQ.seqlen(record1)
        read1::LongDNASeq = FASTQ.sequence(LongDNASeq, record1)
        read1_match_indices::Vector{Int64} = single_match_indices(read1, query_kmer_dict, step)
        sort!(read1_match_indices)
        read1_score::Int64, start_index1, end_index1 = max_window(read1_match_indices, read1_length, length)
        
        if read1_score > threshold
            push!(dataset1_read_numbers, read_num)
            stack_index_range!(coverage, read1_match_indices[start_index1]:read1_match_indices[end_index1])
        end
        
        read2_length = FASTQ.seqlen(record2)
        read2::LongDNASeq = FASTQ.sequence(LongDNASeq, record2)
        read2_match_indices::Vector{Int64} = single_match_indices(read2, query_kmer_dict, step)
        sort!(read2_match_indices)
        read2_score::Int64, start_index2, end_index2 = max_window(read2_match_indices, read2_length, length)
        
        if read2_score > threshold
            push!(dataset2_read_numbers, read_num)
            stack_index_range!(coverage, read2_match_indices[start_index2]:read2_match_indices[end_index2])
        end
    end

    return dataset1_read_numbers, dataset2_read_numbers
end


"Multi-match scanning"
function scan_reads(
    reader1::FASTQ.Reader,
    reader2::FASTQ.Reader,
    query_kmer_dict::Dict{DNAMer{k}, Vector{Int64}},
    read_id_set::Set{Int64},
) where {k}
    nothing
end


function retrieve_reads(
    query_file::String,
    split_dataset::Tuple{String, String};
    output_dir::String="output",
    k::Int64=10,
    threshold::Int64=25,
    step::Int64=1,
    single_match::Bool=false,
)
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    query_sequence, query_length, query_kmer_dict = read_query(query_file, k, true)
    coverage = zeros(Int64, query_length)

    split_dataset_1, split_dataset_2 = split_dataset

    reader1 = FASTQ.Reader(open("$split_dataset_1", "r"))
    reader2 = FASTQ.Reader(open("$split_dataset_2", "r"))

    output_filename = split(split_dataset_1, '_')[1]*".fastq"

    read_numbers_1, read_numbers_2 = scan_reads(reader1, reader2, query_kmer_dict, coverage)

    CSV.write("$output_dir/read_ids.csv", DataFrame(read_ids=sort(collect(union(read_numbers_1, read_numbers_2)))))
    CSV.write("$output_dir/coverage.csv", DataFrame(i=1:query_length, frequency=coverage))
end