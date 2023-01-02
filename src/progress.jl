

struct ProgressBar
    config::Config
    read_count::Int
    file_count::Int
    threads::Int
    progress_meter::ProgressMeter.Progress
end

function ProgressBar(config::Config, read_count::Int, file_count::Int, threads::Int)
    println()
    ProgressBar(
        config, read_count, file_count, threads,
        Progress(read_count, desc="Sniping reads... ", color=:white),
    )
end

function value_list(progress_bar::ProgressBar, reads_scanned::Int, bases_parsed::Int, read_length::Int)
    [
        (Symbol("Reads scanned"), "$reads_scanned/$(progress_bar.read_count)"),
        (Symbol("Bases parsed"), bases_parsed),
        (Symbol("Read length"), read_length),
        (Symbol("k-mer length"), progress_bar.config.k),
        (Symbol("k-mer step"), progress_bar.config.step),
        (Symbol("File count"), progress_bar.file_count),
        (Symbol("Threads utilized"), progress_bar.threads),
    ]
end

function update(progress_bar::ProgressBar, reads_scanned::Int, bases_parsed::Int, read_length::Int, color::Symbol=:white)
    update!(progress_bar.progress_meter, reads_scanned, color, showvalues=value_list(progress_bar, reads_scanned, bases_parsed, read_length))
end

function finish(progress_bar::ProgressBar, bases_parsed::Int, read_length::Int)
    update(progress_bar, progress_bar.read_count, bases_parsed, read_length, :green)
    finish!(progress_bar.progress_meter)
end