
using Plots; gr()


function match_frequency_plot(
    matches::Vector{<:Integer},
    frequency::Vector{<:Integer},
    xlim_upper::Integer = 250,
    step::Integer = 1,
    outfile::AbstractString = "match_frequencies.svg",
)
    N = length(matches)

    plot(matches, (frequency .+ (step-1)) .÷ step,
        color = "#0088ff",
        fillalpha = 0.5,
        fillrange = 1,
        yscale = :log10,
        ylims = (1, Inf),
        xlims = (0, xlim_upper+10),
        xwiden = false,
        tickfont = font(10, "Computer Modern"),
        guidefont = font(12, "Computer Modern"),
        minorgrid = true,
        xlabel = "k-mer matches",
        ylabel = "Frequency",
        legend = false,
        yticks = 8,
        fmt = :svg,
    )

    savefig(outfile)
end


function generate_activity_vectors(reads_SOA, ref_length::Int32, thresholds::Vector{Int64}=[60, 100, 140])

    n = length(reads_SOA)

    matches::Vector{Int32} = reads_SOA.kmer_matches
    range_start::Vector{Int32} = reads_SOA.ref_range_start
    range_end::Vector{Int32} = reads_SOA.ref_range_end

    thresholds = sort(thresholds)
    layer_count = length(thresholds)

    empty_bool_range = fill(false, ref_length)
    activity_vectors::Vector{Vector{Int32}} = [zeros(Int32, ref_length) for _ in 1:layer_count]
    for k in 1:n
        bool_range::Vector{Bool} = map(i -> range_start[k]<=i<=range_end[k], 1:ref_length)
        for (j, (threshold, vec)) in enumerate(zip(thresholds, activity_vectors))
            activity_vectors[j] = vec .+ (threshold < matches[k] ? bool_range : empty_bool_range)
        end
    end
    return activity_vectors
end

export generate_activity_vectors


function activity_plot(activity_vectors::Vector{Vector{<:Integer}}, outfile::AbstractString="activity.png", log::Bool=true)
    
    baseline = zeros(Integer, length(activity_vectors[1]))
    for vec in activity_vectors
        ax.stairs(vec, baseline=baseline, fill=true)
    end
    log ? yscale("log") : yscale("linear")
    savefig(outfile)
end

export activity_plot
