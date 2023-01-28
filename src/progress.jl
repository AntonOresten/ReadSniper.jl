

struct ProgressBar
    config::Config
    file_count::Int
    read_count::Int
    threads::Int
    progress_meter::ProgressMeter.Progress
end

function ProgressBar(config::Config, file_count::Int, read_count::Int, threads::Int)
    ProgressBar(
        config, file_count, read_count, threads,
        Progress(read_count, desc="Sniping reads... ", color=:white),
    )
end

function value_list(progress_bar::ProgressBar, files_scanned::Int, reads_scanned::Int, reads_sniped::Int)
    [
        (Symbol("Files scanned"), "$files_scanned/$(progress_bar.file_count)"),
        (Symbol("Reads scanned"), "$reads_scanned / $(progress_bar.read_count)"),
        (Symbol("Reads sniped"), "$reads_sniped"),
        (Symbol("Threads utilized"), progress_bar.threads),
    ]
end

function update(progress_bar::ProgressBar, files_scanned::Int, reads_scanned::Int, reads_sniped::Int, color::Symbol=:white)
    update!(
        progress_bar.progress_meter,
        reads_scanned,
        color,
        showvalues = value_list(
            progress_bar,
            files_scanned,
            reads_scanned,
            reads_sniped,
        )
    )
end

function finish(progress_bar::ProgressBar, reads_sniped::Int)
    update(progress_bar, progress_bar.file_count, progress_bar.read_count, reads_sniped, :green)
    finish!(progress_bar.progress_meter)
end