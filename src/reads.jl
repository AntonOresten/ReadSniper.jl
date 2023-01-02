

struct ReadResult
    read_number::Int32
    ref_range_start::Int32
    ref_range_end::Int32
    kmer_matches::Int32
end


function highest_score_subseq_in_window(vector::Vector{Int64}, window_size::Int64, score_function::Function=length)
    len = length(vector)
    current_end::Int64 = 0
    best_start::Int64 = 0
    best_end::Int64 = 0
    best_score = 0
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


function scan_read(
    config::Config,
    read_number::Int,
    record::FASTQ.Record,
    query_kmer_dict::Dict{LongSubSeq{DNAAlphabet{4}}, Vector{Int64}},
)
    read_length::Int = seqsize(record)
    read_sequence::LongDNA{4} = sequence(LongDNA{4}, record)
    match_indices_vector::Vector{Vector{Int}} = kmer_match_indices(read_sequence, query_kmer_dict, config.k, config.step)
    
    if length(match_indices_vector) == 0
        return ReadResult(Int32(read_number), zero(Int32), zero(Int32), zero(Int32))
    end
    
    sorted_match_indices::Vector{Int64} = sort(vcat(match_indices_vector...))
    max_matches_in_range::Int64, range_start::Int64, range_end::Int64 = highest_score_subseq_in_window(sorted_match_indices, read_length, length)
    ref_range_start::Int64 = sorted_match_indices[range_start]
    ref_range_end::Int64 = sorted_match_indices[range_end]

    ReadResult(Int32(read_number), Int32(ref_range_start), Int32(ref_range_end), Int32(config.step*max_matches_in_range))
end