

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

    reader = FASTAReader(open(path, "r"))
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

function showinfo(reference::Reference)
    @info """Reference: $(reference.path)
        length: $(reference.length)nt,
        unique k-mers: $(reference.unique_kmer_count)
        GC-content: $(round(100*reference.gc_content, digits=2))%
      """
end


"""
Assumes that the first N reads are representative of the whole Datafile, both in terms of read length and GC content.
Sample size can be set to Inf, such that all read lengths and GC contents get counted.
"""
function get_fasta_metadata(
    path::AbstractString,
    sample_size::Integer = 10000,
    index_occur_count::Integer = 2,
)
    record_counter = 0
    sequence_length_sum = 0
    gc_content_sum = 0.0

    byte_counter = 0
    total_file_size = filesize(path)

    reader = FASTAReader(open(path, "r"))

    for record in reader
        record_counter += 1
        
        byte_counter += record.description_len - (ndigits(record_counter) * index_occur_count) + record.sequence_len + 3 # >\n\n

        sequence_length_sum += record.sequence_len
        gc_content_sum += gc_content(sequence(LongDNA{4}, record))
        if record_counter >= sample_size break end
    end

    close(reader)

    record_size = byte_counter / record_counter

    mean_sequence_length = sequence_length_sum ÷ record_counter
    mean_gc_content = gc_content_sum / record_counter

    if record_counter < sample_size
        record_count = record_counter
    else
        record_count = 0
        for _ in 1:4
            record_count = Int((total_file_size - index_occur_count * cumulative_digit_count(record_count)) ÷ record_size)
        end
    end

    return record_count, mean_sequence_length, mean_gc_content
end


struct DatafileMetadata
    path::AbstractString
    file_size::Int
    read_count::Int
    read_length::Int
    gc_content::Float64
end

function DatafileMetadata(path::AbstractString)::DatafileMetadata
    DatafileMetadata(path, filesize(path), get_fasta_metadata(path)...)
end

function DatafileMetadata(paths::Vector{<:AbstractString})::Vector{DatafileMetadata}
    DatafileMetadata.(paths)
end

export DatafileMetadata

function showinfo(i::Integer, datafile::DatafileMetadata)
    @info """Datafile $i: $(datafile.path)
        size: $(round(datafile.file_size/1e9, digits=3)) GB
        reads: $(datafile.read_count)
        read length: $(datafile.read_length)
        GC-content: $(round(100*datafile.gc_content, digits=2))%
      """
end