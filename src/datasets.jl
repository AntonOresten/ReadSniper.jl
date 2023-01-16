

struct Reference
    path::AbstractString
    sequence::LongDNA{4}
    gc_content::Float64
    length::Int
    kmer_vector::Vector{LongDNA{4}}
    kmer_dict::Dict{LongDNA{4}, Vector{Int}}
    unique_kmer_count::Int
end

function Reference(
    path::AbstractString,
    k::Int,
)::Reference

    reader = FASTA.Reader(open(path, "r"))
    record = first(reader)
    close(reader)

    seq = sequence(LongDNA{4}, record)
    kmer_vector = create_kmer_vector(seq, k)
    kmer_dict = kmer_index_dict(kmer_vector)
    unique_kmer_count = length(unique(keys(kmer_dict)))

    Reference(
        path,
        seq,
        gc_content(seq),
        length(seq),
        kmer_vector,
        kmer_dict,
        unique_kmer_count,
    )
end

export Reference


"""
Assumes that the first N reads are representative of the whole both in terms of read length and GC content.
Sample size can be set to Inf, such that all read lengths and GC contents get counted.
"""
function sample_fasta(path::AbstractString, sample_size::Number=1000)
    read_count::Int = 0

    read_length_sum::Int = 0
    gc_content_sum::Float64 = 0.0

    reader = FASTA.Reader(open(path, "r"))
    
    for record in reader
        read_count += 1

        read_length_sum += seqsize(record)
        gc_content_sum += gc_content(sequence(LongDNA{4}, record))
        if read_count >= sample_size break end
    end
    
    for _ in reader
        read_count += 1
    end

    close(reader)

    mean_read_length = read_length_sum ÷ sample_size
    mean_gc_content = gc_content_sum / sample_size

    return read_count, mean_read_length, mean_gc_content
end


struct FASTAMetadata
    path::AbstractString
    file_size::Int
    read_count::Int
    read_length::Int
    gc_content::Float64
    #index?
end

function FASTAMetadata(path::AbstractString)::FASTAMetadata
    FASTAMetadata(path, filesize(path), sample_fasta(path)...)
end

function FASTAMetadata(paths::Vector{<:AbstractString})::Vector{FASTAMetadata}
    read_count, read_length, gc_content = sample_fasta(paths[1])
    [FASTAMetadata(path, filesize(path), read_count, read_length, gc_content) for path in paths]
end

export FASTAMetadata