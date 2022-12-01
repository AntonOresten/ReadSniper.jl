

using BioSequences, FASTX
using StructArrays

struct ReadMetadata
    id::Int64
    reader::Int64
    match_count::Int64
    range_start::Int64
    range_end::Int64
end


function split_read_sequence(seq::LongDNASeq)
    half_len = length(seq) ÷ 2
    return seq[1:half_len], seq[half_len+1:end]
end


"""
    highest_score_subseq_in_window(
        vector::Vector{Int64},
        window_size::Int64,
        score_function::Function)
"""
function highest_score_subseq_in_window(vector::Vector{Int64}, window_size::Int64, score_function::Function=length)
    len = length(vector)
    current_end::Int64 = 0
    best_start::Int64 = 0
    best_end::Int64 = 0
    best_score::Number = 0
    for current_start::Int64 in 1:len
        while current_end+1 <= len && vector[current_end+1] - vector[current_start] < window_size
            current_end += 1
        end
        current_score = score_function(vector[current_start:current_end])
        if current_score > best_score
            best_start = current_start
            best_end = current_end
            best_score = current_score
        end
    end
    return best_score, best_start, best_end
end


"""
Uses a single-match dictionary when scanning reads.

Works best with small genomes (or large k-values) due to k-mers only being assigned one index.
"""
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

    Threads.@threads for (i, reader) in collect(enumerate((reader1, reader2)))
        for (read_num, record::FASTQ.Record) in ProgressBar(enumerate(reader))
            
            
            if read_score ≥ threshold
                push!(read_ids[i], read_num)
                stack_index_range!(coverage, read_match_indices[start_index]:read_match_indices[end_index])
            end
        end
    end

    return read_ids
end


function analyze_record(
    record::FASTQ.Record,
    query_kmer_dict::Dict{DNAMer{k}, Int64},
    step::Int64,
) where k
    read_length = FASTQ.seqlen(record)
    read::LongDNASeq = FASTQ.sequence(LongDNASeq, record)
    read_match_indices::Vector{Int64} = single_match_indices(read, query_kmer_dict, step)
    read_score::Int64, start_index, end_index = highest_score_subseq_in_window(read_match_indices, read_length, lis)
    

    read_length = FASTQ.seqlen(record)
    read_sequence::LongDNASeq = FASTQ.sequence(LongDNASeq, record)
    read_match_index_vectors::Vector{Vector{Int64}} = multi_match_indices(read_sequence, query_kmer_dict, step)
    concat_indices::Vector{Int64} = vcat(read_match_index_vectors...)
    sorted_indices::Vector{Int64} = sort(concat_indices)
    max_matches_in_range::Int64, range_start::Int64, range_end::Int64 = highest_score_subseq_in_window(sorted_indices, read_length, length)
    return max_matches_in_range, range_start, range_end
end


function analyze_record(
    record::FASTQ.Record,
    query_kmer_dict::Dict{DNAMer{k}, Vector{Int64}},
    step::Int64,
) where k
    read_length = FASTQ.seqlen(record)
    read_sequence::LongDNASeq = FASTQ.sequence(LongDNASeq, record)
    query_index_vectors::Vector{Vector{Int64}} = multi_match_indices(read_sequence, query_kmer_dict, step)
    concat_query_indices::Vector{Int64} = vcat(query_index_vectors...)
    sorted_query_indices::Vector{Int64} = sort(concat_query_indices)
    max_matches_in_range::Int64, range_start::Int64, range_end::Int64 = highest_score_subseq_in_window(sorted_query_indices, read_length, length)
    return max_matches_in_range, range_start, range_end, sorted_query_indices
end
