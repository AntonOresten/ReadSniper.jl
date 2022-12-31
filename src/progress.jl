

struct ProgressBar
    config::Config
    read_count::Int64
    file_count::Int64
    threads::Int64
    progress_meter::ProgressMeter.Progress
end

function ProgressBar( config::Config, read_count::Int64, file_count::Int64, threads::Int64)
    println()
    ProgressBar(
        config, read_count, file_count, threads,
        Progress(read_count, desc="Sniping reads... ", color=:white),
    )
end

function value_list(progress_bar::ProgressBar, reads_scanned::Int64)
    [
        (Symbol("Reads scanned"), "$reads_scanned/$(progress_bar.read_count)"),
        (Symbol("k-mer length"), progress_bar.config.k),
        (Symbol("k-mer step"), progress_bar.config.step),
        (Symbol("File count"), progress_bar.file_count),
        (Symbol("Threads utilized"), progress_bar.threads),
    ]
end

function update(progress_bar::ProgressBar, reads_scanned::Int64, color::Symbol=:white)
    update!(progress_bar.progress_meter, reads_scanned, color, showvalues=value_list(progress_bar, reads_scanned))
end

function finish(progress_bar::ProgressBar)
    update(progress_bar, progress_bar.read_count, :green)
    finish!(progress_bar.progress_meter)
end