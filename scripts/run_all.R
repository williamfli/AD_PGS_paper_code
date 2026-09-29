#!/usr/bin/env Rscript
# scripts/run_all.R -- regenerates every figure panel and table into output_dir (scripts/paths.R)
# and writes manifest.tsv (figure, panel, file, status).
#
#   Rscript scripts/run_all.R                    # full build

CODE_DIR <- local({
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub('^--file=', '', args[grepl('^--file=', args)])
    file_arg <- gsub('~+~', ' ', file_arg, fixed = TRUE)
    script_dir <- if (length(file_arg) == 1) dirname(normalizePath(file_arg)) else getwd()
    normalizePath(file.path(script_dir, '..'))
})
setwd(CODE_DIR)

source('scripts/config.R')
source('src/data.R')
for (module in setdiff(sort(list.files('src', pattern = '\\.R$', full.names = TRUE)), 'src/data.R')) source(module)

pdf(NULL)
QUICK <- nzchar(Sys.getenv('QUICK'))

# Build order; the names are what the command-line filter matches.
BUILD <- list(
    'Associations,TableS3,TableS4,TableS5,TableS6' = make_associations,
    'Figure2'                       = make_figure2,
    'Figure3'                       = make_figure3,
    'Figure4,TableS10'              = make_figure4_prediction,
    'Figure5:pca'                   = make_figure5_pca,
    'Figure5:subtypes'              = make_figure5_subtypes,
    'FigureS1'                      = make_figureS1,
    'FigureS2'                      = make_figureS2,
    'FigureS4'                      = make_figureS4,
    'FigureS5'                      = make_figureS5,
    'FigureS6'                      = make_figureS6,
    'FigureS7'                      = make_figureS7,
    'FigureS8,TableS7'              = make_figureS8,
    'FigureS9,FigureS10'            = make_figureS9_S10,
    'FigureS11'                     = make_figureS11,
    'FigureS12'                     = make_figureS12,
    'FigureS13,TableS11,TableS12'   = if (QUICK) function() make_figureS13_and_tables(n_bootstrap = 50)
                                      else make_figureS13_and_tables,
    'FigureS14'                     = make_figureS14
)

wanted <- commandArgs(trailingOnly = TRUE)
if (length(wanted) > 0) {
    BUILD <- BUILD[vapply(names(BUILD), function(key) any(vapply(wanted, grepl, logical(1), x = key, fixed = TRUE)), logical(1))]
    if (length(BUILD) == 0) stop('run_all.R: no build step matches ', paste(wanted, collapse = ', '))
}

manifest <- list()
failures <- character(0)
for (key in names(BUILD)) {
    say('---- {key}')
    started <- Sys.time()
    rows <- tryCatch(BUILD[[key]](), error = function(e) {
        message('FAILED ', key, ': ', conditionMessage(e))
        failures <<- c(failures, key)
        manifest_row(key, 'ALL', '', 'failed', conditionMessage(e))
    })
    rows$build_step <- key
    rows$seconds <- round(as.numeric(difftime(Sys.time(), started, units = 'secs')), 1)
    manifest[[key]] <- rows
    say('     {nrow(rows)} output(s) in {rows$seconds[1]}s')
}
manifest <- do.call(rbind, manifest)
manifest$path <- sub(OUTPUT_DIR, '', manifest$path, fixed = TRUE)
manifest$built_at <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')
manifest$r_version <- R.version.string

manifest_file <- glue('{OUTPUT_DIR}{if (length(wanted) > 0) "manifest_partial" else "manifest"}.tsv')
write.table(manifest, manifest_file, sep = '\t', quote = FALSE, row.names = FALSE)

cat('\n==== Build summary ====\n')
print(table(manifest$status))
cat('total time:', round(sum(manifest$seconds[!duplicated(manifest$build_step)]) / 60, 1), 'min\n')
cat('manifest:  ', manifest_file, '\n')
invisible(dev.off())
if (length(failures) > 0) {
    cat('FAILED steps: ', paste(failures, collapse = ', '), '\n')
    quit(status = 1)
}
