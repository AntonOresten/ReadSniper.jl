
function csv_to_read_index_set(
    csv_file::AbstractString,
)
    df = CSV.read(csv_file, DataFrame)
    Set(df.read_index)
end


function filter_fasta(
    read_index_set::Set{<:Integer},
    fasta_paths::Vector{<:AbstractString},
    output_file::AbstractString,
)
    readers = FASTAReader.(open.(fasta_paths, "r"))
    writer = FASTAWriter(open(output_file, "w"))

    for (read_index, records) in enumerate(zip(readers...))
        if read_index in read_index_set
            foreach(record->write(writer, record), records)
        end
    end

    close.(readers)
    close(writer)
end


export read_index_set, filter_fasta