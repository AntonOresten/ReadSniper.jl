

struct ProgressBar
    config::Config
    file_count::Int
    read_count::Int
    progress_meter::ProgressMeter.Progress
end

function ProgressBar(config::Config, file_count::Int, read_count::Int)
    ProgressBar(
        config, file_count, read_count,
        Progress(read_count, desc="Sniping reads... ", color=:white),
    )
end

function value_list(progress_bar::ProgressBar, files_scanned::Int, reads_scanned::Int, read_count::Int, reads_sniped::Int)
    [
        (Symbol("Files scanned"), "$files_scanned/$(progress_bar.file_count)"),
        (Symbol("Reads scanned"), "$reads_scanned / $read_count"),
        (Symbol("Reads sniped"), "$reads_sniped"),
    ]
end

function update(
    progress_bar::ProgressBar,
    files_scanned::Int,
    reads_scanned::Int,
    reads_sniped::Int,
    color::Symbol=:white;
    finished::Bool=false,
)
    update!(
        progress_bar.progress_meter,
        finished ? progress_bar.read_count : min(reads_scanned, progress_bar.read_count-1),
        color,
        showvalues = value_list(
            progress_bar,
            files_scanned,
            reads_scanned,
            finished ? reads_scanned : progress_bar.read_count,
            reads_sniped,
        )
    )
end

function finish(progress_bar::ProgressBar, reads_scanned::Int, reads_sniped::Int)
    update(progress_bar, progress_bar.file_count, reads_scanned, reads_sniped, :green, finished=true)
    finish!(progress_bar.progress_meter)
end