using Pkg
#Pkg.add(PackageSpec(name="ReadSniper", url="https://github.com/Periareion/ReadSniper.jl"))

using ReadSniper
using Test

@testset "ReadSniper.jl" begin
    
    @test lis([1,7,14,4,5,4,6,15]) == 5
    @test lis([1,7,14,4,5,5,6,15]) == 5
end
