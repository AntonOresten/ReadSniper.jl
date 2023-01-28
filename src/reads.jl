

struct ReadResult
    read_number::Int32
    ref_range_start::Int32
    ref_range_end::Int32
    score::Int32
end


function scan_read(
    config::Config,
    record::FASTARecord,
    kmer_dict::Dict{LongDNA{4}, Vector{Int64}},
)::Tuple{Int,Int,Int}

    read_length::Int = seqsize(record)
    read_sequence::LongDNA{4} = sequence(LongDNA{4}, record)
    match_indices::Vector{Int} = vcat(kmer_match_indices(read_sequence, kmer_dict, config.k, config.step)...)

    if length(match_indices) == 0
        return (0, 0, 0)
    end

    sort!(match_indices)
    score::Int64, range_start::Int64, range_end::Int64 = max_subseq_in_range(match_indices, read_length)
    ref_range_start::Int64 = match_indices[range_start]
    ref_range_end::Int64 = match_indices[range_end]

    if ref_range_start == -769 && ref_range_end == -545
        scatter(
            match_indices,
            1:length(match_indices),
            tickfont = font(10, "Computer Modern"),
            guidefont = font(12, "Computer Modern"),
            xlabel="Position in reference",
            formatter=:plain,
            legend = false,
            fmt = :svg,
            margins = 0.4Plots.cm,
        )
        savefig("read.svg")
    end

    return (ref_range_start * !(ref_range_start < 0 < ref_range_end), ref_range_end, config.step*score)
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
    #     push!(read_results, ReadResult(Int32(read_number), Int32(mw[1]), Int32(mw[2]), Int32(config.step*score)))
    # end

    # read_results