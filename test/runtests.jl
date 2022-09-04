using Pkg
Pkg.add(PackageSpec(name="ReadSniper", url="https://github.com/Periareion/ReadSniper.jl"))


using ReadSniper
using Test

@show max_window([1, 2, 3, 4, 5, 6, 100, 110], 5)
@show max_window([1, 2, 5, 6, 100, 110, 109, 108, 107, 106, 105, 106, 107, 108, 109, 110], 5)

@testset "ReadSniper.jl" begin
    # Write your tests here.
end
