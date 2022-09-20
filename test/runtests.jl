using Pkg
#Pkg.add(PackageSpec(name="ReadSniper", url="https://github.com/Periareion/ReadSniper.jl"))

include("c:\\Users\\anton\\Code\\Local-Repositories\\ReadSniper.jl\\src\\ReadSniper.jl")
using Test

@testset "ReadSniper.jl" begin
    
    @test ReadSniper.lis([1,7,14,4,5,4,6,15]) == 5
    @test ReadSniper.lis([1,7,14,4,5,5,6,15]) == 5
end