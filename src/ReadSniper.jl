
module ReadSniper

    using Dates
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
    include("filter.jl")
    include("probability.jl")
    include("iteration.jl")
    include("plots.jl")
    include("main.jl")

end

# Hello, person reading this in the future (probably me)! :)
# Don't forget to add capability to snipe reads in FASTQ files.