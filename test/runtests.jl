
using ReadSniper
using Test

dir = "C:/Users/anton/RSData"

retrieve_reads(
    dir*"/queries/picorna_nuevo_rc.fasta",
    (
        dir*"/datasets/SRR10873757/SRR10873757_1.fastq",
        dir*"/datasets/SRR10873757/SRR10873757_2.fastq"
    ),
    output_dir="output",
    k=10,
    step=1,
    single_match=false
)


#reads_fields_tuple = CSV.read("output2/reads.csv", NamedTuple)
#reads_fields_tuple = (max_matches_in_range=[15,25,35,25], range_start=[10,40,70,20], range_end=[60, 80, 90, 60])
#activity_vectors = generate_activity_vectors(reads_fields_tuple, 20000)
#println(activity_vectors)
#plot_activity_histogram(activity_vectors)


#=
using CSV

df = CSV.read("output/read_match_distribution.csv", NamedTuple)

scores = df.scores
counts = df.counts

bar_plot(scores, counts, "read_match_distribution.png", true)
=#


#=
@testset "ReadSniper.jl" begin
    @test lis([1,7,14,4,5,4,6,15]) == 5
    #@test lis([20,1,7,14,4,5,5,6,15]) == 5
    #@test single_match_indices(dna"ACGGTCGATCTG", single_match_kmer_dict(dna"ACGAGTCGATTCCTG", 3)) == [1, 5, 6, 2, 8, 13]
    #single_match_kmer_dict(dna"ACGAGTCGATTCCTG", 3)
end
=#
