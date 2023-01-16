
using ReadSniper
using Test

using BioSequences


@testset "utils.jl" begin
    @test lis([1,7,14,4,5,4,6,15]) == 5
    @test lis([20,1,7,14,4,5,5,6,15]) == 5
    
    dict = Dict(1 => 1, 2 => 1)
    increment_dict_value!(1, dict)
    @test dict == Dict(1 => 2, 2 => 1)
end


@testset "probability.jl" begin
    nothing
end


@testset "kmers.jl" begin
    @test create_kmer_vector(dna"ACGTACGT", 4) == [dna"ACGT", dna"CGTA", dna"GTAC", dna"TACG", dna"ACGT"]

    @testset "kmer index dictionaries" begin
        @test kmer_index_dict(dna"ACGTACGT", 4, false) == Dict(dna"ACGT" => [5, 1], dna"CGTA" => [2], dna"GTAC" => [3], dna"TACG" => [4])
        @test kmer_match_indices(dna"ACGTTCGTA", kmer_index_dict(dna"ACGTACGT", 4, false), 4) == [[5, 1], [2]]
    end
end


@testset "reads.jl" begin
    @test highest_score_subseq_in_window([1,6,7,8,9,10,11,12,13,17], 8, length) == (8, 2, 9)
end


@testset "datasets.jl" begin
    # reference
    # fastqfile
    nothing
end


if isdir("C:/Users/anton/RSData")
    snipe_reads(
        "C:/Users/anton/RSData/reference/reference.fasta",
        ["SRR10873757_1.fasta", "SRR10873757_2.fasta"],
        "C:/Users/anton/RSData/datasets/fasta-files/SRR10873757",
        output_dir="output",
        k=9, step=4,
    )
end
