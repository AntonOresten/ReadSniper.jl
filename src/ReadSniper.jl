
module ReadSniper

    using BioSequences, FASTX
    using CSV, StructArrays
    using LinearAlgebra
    using ProgressMeter

    include("utils.jl")
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

# write all read scores to a file maybe - 16 bits per score, and then enumerate when reading to get read num