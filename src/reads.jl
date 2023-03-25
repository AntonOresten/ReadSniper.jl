

struct ReadResult
    read_index::Int32
    ref_range_start::Int32
    ref_range_end::Int32
    score::Int32
end


function scan_read(
    config::Config,
    record::FASTARecord,
    kmer_dict::Dict{LongDNA{4}, Vector{Int}},
    read_num::Int,
)::Tuple{Int,Int,Int}

    read_length::Int = seqsize(record)
    read_sequence::LongDNA{4} = sequence(LongDNA{4}, record)
    index_vectors::Vector{Vector{Int}} = kmer_match_indices(read_sequence, kmer_dict, config.k, config.step)

    if length(index_vectors) == 0 return (0, 0, 0) end

    score::Int, ref_range_start::Int, ref_range_end::Int = constrained_multi_choice_LIS(index_vectors, read_length)
    
    return score*config.step, ref_range_start * !(ref_range_start < 0 < ref_range_end), ref_range_end
end

export scan_read


    # for match_indices in match_indices_vector
    #     for index in match_indices
    #         bool_match_vector[abs(index)] = true
    #     end
    # end

    # matched_windows = remove_overlapping_pairs(find_subarrays(bool_match_vector, read_length-config.k+1, config.threshold ÷ config.step))
    
    # read_results::Vector{ReadResult} = []

    # for mw in matched_windows
    #     score = count(view(bool_match_vector, mw[1]:mw[2]))
    #     push!(read_results, ReadResult(Int32(read_index), Int32(mw[1]), Int32(mw[2]), Int32(config.step*score)))
    # end

    # read_results