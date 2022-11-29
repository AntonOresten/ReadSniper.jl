
using ReadSniper
using Test

dir = "C:/Users/anton/RSData"

retrieve_reads(
    dir*"/queries/picorna_nuevo_rc.fasta",
    (
        dir*"/datasets/SRR10873757/SRR10873757_1.fastq",
        dir*"/datasets/SRR10873757/SRR10873757_2.fastq"
    ),
    k=10,
    threshold=25,
    step=1,
    single_match=false
)

#=
@testset "ReadSniper.jl" begin
    @test lis([1,7,14,4,5,4,6,15]) == 5
    #@test lis([20,1,7,14,4,5,5,6,15]) == 5
    #@test single_match_indices(dna"ACGGTCGATCTG", single_match_kmer_dict(dna"ACGAGTCGATTCCTG", 3)) == [1, 5, 6, 2, 8, 13]
    #single_match_kmer_dict(dna"ACGAGTCGATTCCTG", 3)
end
=#