
using PyPlot


function bar_plot(x::Vector, y::Vector, outfile::String, log::Bool=false)
    plt.style.use("seaborn")
    p = bar(x, y, width=1)
    if log
        yscale("log")
    end
    savefig(outfile, dpi=1000)
end


function score_frequency_plot(scores::Vector{Int32}, counts::Vector{Int32}, outfile::String="score_frequencies.png")
    bar_plot(scores, counts, outfile, true)
end

export score_frequency_plot


function generate_activity_vectors(reads_SOA, query_length::Int32, thresholds::Vector{Int64}=[60, 100, 140])
    n = length(reads_SOA.ref_range_start)

    score::Vector{Int32} = reads_SOA.kmer_match_count
    range_start::Vector{Int32} = reads_SOA.ref_range_start
    range_end::Vector{Int32} = reads_SOA.ref_range_end

    thresholds = sort(thresholds)
    layer_count = length(thresholds)

    activity_vectors::Vector{Vector{Int32}} = [zeros(Int32, query_length) for _ in 1:layer_count]
    for k in 1:n
        binary_range = map(i -> range_start[k]<=i<=range_end[k], 1:query_length)
        for (j, (threshold, vec)) in enumerate(zip(thresholds, activity_vectors))
            if threshold > score[k]
                break
            end
            activity_vectors[j] = vec .+ binary_range
        end
    end
    return activity_vectors
end

export generate_activity_vectors


function activity_plot(activity_vectors::Vector{Vector{Int32}}, outfile::String="activity.png", log::Bool=true)
    plt.style.use("seaborn")

    fig, ax = subplots()
    
    baseline = zeros(Int32, length(activity_vectors[1]))
    for vec in activity_vectors
        ax.stairs(vec, baseline=baseline, fill=true)
    end
    log ? yscale("log") : nothing
    savefig(outfile, dpi=1000)
end

export activity_plot
