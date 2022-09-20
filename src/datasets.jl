
using CSV
using DataFrames

function get_best_datasets(sra_metadata_file::String, n::Int64)
    datasets::DataFrame = DataFrames.sort(CSV.read(sra_metadata_file, DataFrame), [:pident], rev=true)
    run_ids = datasets.run_id
    best_SRR::Vector{String} = []
    for run_id in run_ids
        if length(best_SRR) >= n
            break
        end
        if occursin("SRR", run_id)
            push!(best_SRR, run_id)
        end
    end
    return best_SRR[1:n]
end

function download_datasets(datasets::Vector{String})
    nothing
end

@show get_best_datasets("src/_tmp_submission.sra.metadata.csv", 20)