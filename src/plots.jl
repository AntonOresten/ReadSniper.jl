
using Plots; gr()
using Colors


function init_score_distribution_plot(
    xlim_upper::Integer = 250,
)
    plot(
        yscale = :log10,
        yticks = 10,
        ylims = (1, Inf),
        xlims = (0, xlim_upper+10),
        tickfont = font(10, "Computer Modern"),
        guidefont = font(12, "Computer Modern"),
        minorgrid = true,
        xlabel = "Read score",
        ylabel = "Frequency",
        yminorgrid = true,
        legendfont = font(10, "Computer Modern"),
        legend_position = :topright, 
        fmt = :svg,
        #margins = 0.5Plots.cm,
    )
end

function score_distribution_plot(
    matches::Vector{<:Integer},
    frequency::Vector{<:Integer},
    k::Integer,
    step::Integer,
    xlim_upper::Integer = 250,
    estimate::Union{Vector{<:AbstractFloat}, Nothing} = [],
    outfile::AbstractString = "score_distribution.svg";
    result_color = "#0088ff",
    estimate_color = "#ff8800",
)
    name = "k$(k)s$(step)"

    init_score_distribution_plot(xlim_upper)

    adjusted_frequency = (frequency .+ (step-1)) .÷ step

    plot!(matches, adjusted_frequency,
        color = result_color,
        fillalpha = 0.5,
        fillrange = 1,
        label = "Empirical ($name)",
    )

    if estimate !== nothing
        values = map(x->max(1, x), estimate)
        plot!(
            filter(p -> p[2] > 1, collect(enumerate(values))),
            color = estimate_color,
            fillalpha = 0.5,
            fillrange = 1,
            label = "Expected ($name)"
        )
    end

    savefig(outfile)
end

function score_distribution_plot(
    file::AbstractString,
    k::Integer,
    step::Integer = 1,
    xlim_upper::Integer = 250,
    estimate::Union{Vector{<:AbstractFloat}, Nothing} = [],
    outfile::AbstractString = "score_distribution.svg",
)
    match_frequencies = CSV.read(file, NamedTuple)
    score_distribution_plot(
        match_frequencies.matches,
        match_frequencies.frequency,
        k,
        step,
        xlim_upper,
        estimate,
        outfile,
    )
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
    
    log_mode_min = 0.9
    plot(
        tickfont = font(10, "Computer Modern"),
        guidefont = font(12, "Computer Modern"),
        xlabel = "Position in reference",
        ylabel = "Match activity",
        yscale = log_scale ? :log10 : :identity,
        ylim = log_scale ? (log_mode_min, Inf) : (0, Inf),
        yminorgrid = true,
        legendfont = font(10, "Computer Modern"),
        fmt = :svg,
        legend = :outertopright,
        xformatter = :plain,
        xticks = [1, 2500, 5000, 7500, 9366],
    )
    for (i, vec) in enumerate(activity_vectors)
        values = map(x->max(0.1, x), vec)
        plot!(
            values,
            color = colors[i],
            fillalpha = 1,
            fillrange = log_scale ? log_mode_min : 0,
            yscale = log_scale ? :log10 : :identity,
            label = "$(thresholds[i])",
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