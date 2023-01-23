
using Plots; gr()
using Colors


function init_score_distribution_plot(
    xlim_upper::Integer = 250,
)
    plot(
        yscale = :log10,
        ylims = (1, Inf),
        xlims = (0, xlim_upper+10),
        tickfont = font(10, "Computer Modern"),
        guidefont = font(12, "Computer Modern"),
        minorgrid = true,
        xlabel = "Read score",
        ylabel = "Frequency",
        yticks = 4,
        fmt = :svg,
    )
end

function score_distribution_plot(
    matches::Vector{<:Integer},
    frequency::Vector{<:Integer},
    k::Integer,
    step::Integer,
    estimate::Union{Vector{<:AbstractFloat}, Nothing} = [],
    color1 = "#0088ff",
    color2 = "#ff8800",
)
    name = "k$(k)s$(step)"

    plot!(matches, (frequency .+ (step-1)) .÷ step,
        color = color1,
        fillalpha = 0.5,
        fillrange = 1,
        label = "Result ($name)",
    )

    if estimate !== nothing
        plot!(
            map(x->max(1, x), estimate),
            color = color2,
            fillalpha = 0.5,
            fillrange = 1,
            label = "Estimate ($name)"
        )
    end
end

function score_distribution_plot(
    file::AbstractString,
    k::Integer,
    step::Integer = 1,
    xlim_upper::Integer = 250,
    estimate::Union{Vector{<:AbstractFloat}, Nothing} = [],
    outfile::AbstractString = "score_distribution.svg",
)
    init_score_distribution_plot(xlim_upper)

    match_frequencies = CSV.read(file, NamedTuple)
    score_distribution_plot(match_frequencies.matches, match_frequencies.frequency, k, step, estimate)

    savefig(outfile)
end

export score_distribution_plot


function generate_activity_vectors(read_result_fields, ref_length::Integer, k::Integer, thresholds::Vector{Int64}=[60, 100, 140])

    matches::Vector{Int32} = read_result_fields.score
    range_start::Vector{Int32} = read_result_fields.ref_range_start
    range_end::Vector{Int32} = read_result_fields.ref_range_end

    n = length(matches)

    thresholds = sort(thresholds)
    layer_count = length(thresholds)
    read_counts = zeros(Int, layer_count)

    activity_vectors::Vector{Vector{Int32}} = [zeros(Int32, ref_length) for _ in 1:layer_count]
    for i in 1:n
        A, B = range_start[i] > 0 ? (range_start[i], range_end[i]+k-1) : (-range_end[i], -range_start[i]+k-1)
        bool_range::Vector{Bool} = map(x -> A<=x<=B, 1:ref_length)
        for (j, threshold) in enumerate(thresholds)
            if matches[i] >= threshold
                read_counts[j] += 1
                activity_vectors[j] += bool_range
            end
        end
    end
    return activity_vectors
end

export generate_activity_vectors


function activity_plot(
    activity_vectors::Vector{<:Vector{<:Integer}},
    thresholds::Vector{<:Integer},
    outfile::AbstractString="activity.svg",
    log_scale::Bool=false,
)
    N = length(activity_vectors)
    colors = range(colorant"plum", colorant"darkmagenta", length=N)

    max_activity = maximum(maximum.(activity_vectors))
    
    plot(
        tickfont = font(10, "Computer Modern"),
        guidefont = font(12, "Computer Modern"),
        xlabel = "Position in reference",
        ylabel = "Match activity",
        yscale = log_scale ? :log10 : :identity,
        ylims = log_scale ? (1, max_activity) : (0, max_activity),
        yticks = 3,
        fmt = :svg,
        legend  = :outertopright
    )
    for (i, vec) in enumerate(activity_vectors)
        plot!(
            color = colors[i],
            map(x->max(1, x), vec),
            fillalpha = 1,
            fillrange = log_scale ? 1 : 0,
            yscale = log_scale ? :log10 : :identity,
            label = ">$(thresholds[i]) score",
        )
    end
    savefig(outfile)
end

function activity_plot(
    read_result_fields,
    ref_length::Integer,
    k::Integer,
    thresholds::Vector{<:Integer},
    outfile::AbstractString,
    log_scale::Bool=false,
)
    activity_plot(
        generate_activity_vectors(
            read_result_fields,
            ref_length,
            k,
            thresholds,
        ),
        thresholds,
        outfile,
        log_scale
    )
end

export activity_plot