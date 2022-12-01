
module ReadSniper

    using BioSequences, FASTX
    using CSV, StructArrays
    using ProgressMeter
    import Dates

    include("utils.jl")
    include("plots.jl")
    include("kmers.jl")
    include("datasets.jl")
    include("read_iteration.jl")
    include("reads.jl")

end