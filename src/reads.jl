

struct ReadMetadata
    id::Int64
    reader::Int64
    match_count::Int64
    range_start::Int64
    range_end::Int64
end


struct DatafileReader
    reader::FASTQ.Reader
    id::Int64
    window::Int64
    first_index::Int64
    last_index::Int64
end

calc_first_index(id, window, remainder) = window * (id-1) + 1
calc_last_index(id, window, remainder) = window * id # + (id == 1 ? remainder : 0)

function DatafileReader(reader::FASTQ.Reader, id::Int64, window::Int64, remainder::Int64)
    first_index = calc_first_index(id, window, remainder)
    last_index = calc_last_index(id, window, remainder)
    df_reader = DatafileReader(reader, id, last_index-first_index+1, first_index, last_index)
    return skip_to_window(df_reader)
end


function skip_to_window(datafile_reader::DatafileReader)
    reader = datafile_reader.reader
    for _ in 1:datafile_reader.first_index-1
        iterate(reader)
    end
    return datafile_reader
end


function split_read_sequence(seq::LongDNA{4})
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

export highest_score_subseq_in_window


"""
Uses a single-match dictionary when scanning read

Works best with small genomes (or large k-values) due to k-mers only being assigned one index.
"""
function analyze_record(
    record::FASTQ.Record,
    query_kmer_dict::Dict{LongSubSeq{DNAAlphabet{4}}, Int64},
    k_1::Int64,
    step::Int64,
    threshold::Int64,
)
    read_length::Int64 = seqsize(record)
    read_sequence::LongDNA{4} = FASTQ.sequence(LongDNA{4}, record)
    read_match_indices::Vector{Int64} = single_match_indices(read_sequence, query_kmer_dict, k_1, step)
    sort!(read_match_indices)
    if length(read_match_indices) < threshold
        return 0, 1, 1, [1]
    else
        max_matches_in_range::Int64, range_start::Int64, range_end::Int64 = highest_score_subseq_in_window(read_match_indices, read_length, length)
        return max_matches_in_range, range_start, range_end, read_match_indices
    end
end

"""
Uses a multi-match dictionary when scanning read
"""
function analyze_record(
    record::FASTQ.Record,
    query_kmer_dict::Dict{LongSubSeq{DNAAlphabet{4}}, Vector{Int64}},
    k_1::Int64,
    step::Int64,
    threshold::Int64,
)
    read_length::Int64 = seqsize(record)
    read_sequence::LongDNA{4} = FASTQ.sequence(LongDNA{4}, record)
    read_match_index_vectors::Vector{Vector{Int64}} = multi_match_indices(read_sequence, query_kmer_dict, k_1, step)
    
    if length(read_match_index_vectors) < threshold
        return 0, 1, 1, [1]
    else
        sorted_concat_indices::Vector{Int64} = sort(vcat(read_match_index_vectors...))
        max_matches_in_range::Int64, range_start::Int64, range_end::Int64 = highest_score_subseq_in_window(sorted_concat_indices, read_length, length)
        return max_matches_in_range, range_start, range_end, sorted_concat_indices
    end
end
