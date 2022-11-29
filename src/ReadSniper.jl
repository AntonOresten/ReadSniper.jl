
module ReadSniper

    using BioSequences, FASTX
    using CSV, DataFrames
    using ProgressBars
    import Dates

    include("utils.jl")
    #include("plots.jl")
    include("kmers.jl")
    include("datasets.jl")
    include("reads.jl")

end