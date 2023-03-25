

"""
Naive approximation (assumes larger k and not too large L).
A better approximation would account for duplicate kmers and GC content.
"""
kmer_count(L::Int, k::Int)::Int = L - k + 1

nucleotide_combinations(k::Int)::Int = 4^k

P_ACGT(GC::AbstractFloat)::Vector{<:AbstractFloat} = [1-GC, GC, GC, 1-GC] .* 0.5

P_nucleotide_match(GC1::AbstractFloat, GC2::AbstractFloat)::AbstractFloat = P_ACGT(GC1) ⋅ P_ACGT(GC2)

P_kmer_match(k::Integer, GC1::AbstractFloat, GC2::AbstractFloat)::AbstractFloat = P_nucleotide_match(GC1, GC2)^k


function binomial_probability(n::Int, x::Int, p::AbstractFloat)
    binomial(BigInt(n), BigInt(x))*p^x*(1-p)^(n-x)
end


function estimate_score_distribution(
    k::Int,
    reference::Reference,
    datafile_list::Vector{DatafileMetadata},
)::Vector{<:AbstractFloat}
    total_read_count = sum([fqm.read_count for fqm in datafile_list])
    datasets_gc_content = sum([fqm.gc_content for fqm in datafile_list]) / length(datafile_list)
    read_length = datafile_list[1].read_length

    # The important thing here is to have an approximation of a random score distribution that scales correctly with k and read length.
    p = 1 - (1 - P_kmer_match(k, reference.gc_content, datasets_gc_content)*sqrt(read_length)/150*(k > 5 ? (k - 5)^1.4 : 1))^reference.unique_kmer_count
    #total_read_count * binomial_probability.(read_length, 0:read_length, p)

    #p = 1 - (1 - (1 - (1 - P_kmer_match(k, reference.gc_content, datasets_gc_content)*(read_length-k+1)/(reference.length - k + 1))^(2*(reference.length - k + 1)/reference.unique_kmer_count)))^(reference.unique_kmer_count)
    # * 2 * (reference.length - k + 1) ÷ reference.unique_kmer_count
    total_read_count * binomial_probability.(read_length, 0:read_length, p)
end

export estimate_score_distribution

#println(findfirst(map(x->x[2]<0&&x[1]>15,enumerate(log10.(estimate_match_distribution(8, 9360, 246, 2575656))))))

export P_ACGT, P_kmer_match