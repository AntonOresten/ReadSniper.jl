

"""
Naive approximation (assumes larger k and not too large L).
A better approximation would account for duplicate kmers and GC content.
"""
kmer_count(L::Int, k::Int)::Int = L - k + 1

nucleotide_combinations(k::Int)::Int = 4^k

P_ACGT(GC::AbstractFloat)::Vector{<:AbstractFloat} = [1-GC, GC, GC, 1-GC] / 2

P_nucleotide_match(GC1::AbstractFloat, GC2::AbstractFloat)::AbstractFloat = P_ACGT(GC1) ⋅ P_ACGT(GC2)

P_kmer_match(k::Integer, GC1::AbstractFloat, GC2::AbstractFloat)::AbstractFloat = P_nucleotide_match(GC1, GC2)^k


function binomial_probability(n::Int, x::Int, p::AbstractFloat)
    binomial(BigInt(n), BigInt(x))*p^x*(1-p)^(n-x)
end


function spurious_score_distribution(
    k::Int,
    reference::Reference,
    fasta_metadata_list::Vector{FASTAMetadata},
)::Vector{<:AbstractFloat}
    total_read_count = sum([fqm.read_count for fqm in fasta_metadata_list])
    datasets_gc_content = sum([fqm.gc_content for fqm in fasta_metadata_list]) / length(fasta_metadata_list)
    read_length = fasta_metadata_list[1].read_length

    # The important thing here is to have an approximation of a random score distribution that scales correctly with k and read length.
    p = 1 - (1 - P_kmer_match(k, reference.gc_content, datasets_gc_content)/7)^reference.unique_kmer_count
    #total_read_count * binomial_probability.(read_length, 0:read_length, p)

    #p = 1 - (1 - (1 - (1 - P_kmer_match(k, reference.gc_content, datasets_gc_content)*(read_length-k+1)/(reference.length - k + 1))^(2*(reference.length - k + 1)/reference.unique_kmer_count)))^(reference.unique_kmer_count)
    HIGH = read_length# * 2 * (reference.length - k + 1) ÷ reference.unique_kmer_count
    total_read_count * binomial_probability.(HIGH, 0:HIGH, p)
end


#println(findfirst(map(x->x[2]<0&&x[1]>15,enumerate(log10.(spurious_match_distribution(8, 9360, 246, 2575656))))))
