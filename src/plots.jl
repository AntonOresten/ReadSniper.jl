
using PyPlot


function bar_plot(x::Vector, y::Vector, save_as::String="", log::Bool=false)
    p = bar(x, y, width=1)
    if log
        yscale("log")
    end
    savefig(save_as, dpi=1000)
end


function generate_activity_vectors(reads_fields_tuple::NamedTuple, query_length::Int64, thresholds::Vector{Int64}=[10,50,100])
    n = length(reads_fields_tuple.range_start)

    max_matches_in_range::Vector{Int64} = reads_fields_tuple.max_matches_in_range
    range_start::Vector{Int64} = reads_fields_tuple.range_start
    range_end::Vector{Int64} = reads_fields_tuple.range_end

    thresholds = sort(thresholds)
    layer_count = length(thresholds)

    activity_vectors::Vector{Vector{Int64}} = [zeros(Int64, query_length) for _ in 1:layer_count]
    for k in 1:n
        binary_range = map(i -> range_start[k]<=i<=range_end[k], 1:query_length)
        #println(collect(binary_range))
        for (j, (threshold, vec)) in enumerate(zip(thresholds, activity_vectors))
            if threshold > max_matches_in_range[k]
                break
            end
            activity_vectors[j] = vec .+ binary_range
            #println("addind $vec")
        end
    end
    return activity_vectors
end

export generate_activity_vectors


function plot_activity_histogram(activity_vectors::Vector{Vector{Int64}}, save_as::String="activity.png")
    plt.style.use("seaborn")

    fig, ax = subplots()
    
    baseline = zeros(Int64, length(activity_vectors[1]))
    for vec in activity_vectors
        println("meh")
        ax.stairs(vec, baseline=baseline, fill=true)
    end
    show()
    #savefig(save_as, dpi=1000)
end

export plot_activity_histogram

#3plot_coverage_histogram()


#=
df = CSV.read("test/output/coverage.csv", DataFrame)

x, y = df.i, df.frequency
n = length(x)

bar(x, y, width=1)
show()
=#