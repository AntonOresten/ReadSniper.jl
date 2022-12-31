
module ReadSniper

    using BioSequences, FASTX
    using CSV, StructArrays
    using ProgressMeter
    import Dates

    include("utils.jl")
    include("config.jl")
    include("plots.jl")
    include("kmers.jl")
    include("reads.jl")
    include("datasets.jl")
    include("progress.jl")
    include("read_iteration.jl")
    include("scoring.jl")
    include("API.jl")
end
