
module ReadSniper

    using BioSequences, FASTX
    using CSV, StructArrays
    using LinearAlgebra
    using ProgressMeter

    include("utils.jl")
    include("randutils.jl")
    include("config.jl")
    include("kmers.jl")
    include("reads.jl")
    include("progress.jl")
    include("datasets.jl")
    include("probability.jl")
    include("iteration.jl")
    include("plots.jl")
    include("main.jl")

end