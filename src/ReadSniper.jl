
module ReadSniper

    using BioSequences, FASTX
    using CSV, StructArrays
    using ProgressMeter

    include("utils.jl")
    include("config.jl")
    include("plots.jl")
    include("kmers.jl")
    include("reads.jl")
    include("progress.jl")
    include("iteration.jl")
    include("datasets.jl")
    include("API.jl")

end
