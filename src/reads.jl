

struct ReadResult
    read_number::Int32
    ref_range_start::Int32
    ref_range_end::Int32
    kmer_matches::Int32
end


function scan_read(
    config::Config,
    read_number::Int,
    record::FASTARecord,
    kmer_dict::Dict{LongDNA{4}, Vector{Int64}},
    #bool_match_vector::Vector{Bool},
)
    read_length::Int = seqsize(record)
    read_sequence::LongDNA{4} = sequence(LongDNA{4}, record)
    match_indices_vector::Vector{Vector{Int}} = kmer_match_indices(read_sequence, kmer_dict, config.k, config.step)
    
    if length(match_indices_vector) == 0
        return ReadResult(Int32(read_number), Int32(0), Int32(0), Int32(0))
    end
    
    sorted_match_indices::Vector{Int64} = sort(vcat(match_indices_vector...))
    score::Int64, range_start::Int64, range_end::Int64 = highest_score_subseq_in_window(sorted_match_indices, read_length, length)
    ref_range_start::Int64 = sorted_match_indices[range_start]
    ref_range_end::Int64 = sorted_match_indices[range_end]

    ReadResult(Int32(read_number), Int32(ref_range_start), Int32(ref_range_end), Int32(config.step*score))

    # for match_indices in match_indices_vector
    #     for index in match_indices
    #         bool_match_vector[abs(index)] = true
    #     end
    # end

    # matched_windows = remove_overlapping_pairs(find_subarrays(bool_match_vector, read_length-config.k+1, config.threshold))
    
    # read_results::Vector{ReadResult} = []

    # for (i, mw) in enumerate(matched_windows)
    #     score = count(view(bool_match_vector, mw[1]:mw[2]))
    #     push!(read_results, ReadResult(Int32(read_number), Int32(mw[1]), Int32(mw[2]), Int32(config.step*score)))
    #     if i == 1
    #         break
    #     end
    # end

    # read_results
    

end