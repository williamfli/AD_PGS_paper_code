# scripts/config.R -- shared settings for the figure build (analysis constants and helpers). Sourced by scripts/run_all.R.
suppressPackageStartupMessages({
    library(glue)
    library(ggplot2)
})

# ---- Paths -----------------------------------------------------------------------------
PATHS_FILE <- Sys.getenv('AD_PGS_PATHS_FILE')
if (!nzchar(PATHS_FILE)) PATHS_FILE <- 'scripts/paths.R'
if (!file.exists(PATHS_FILE)) stop(glue('config.R: paths file not found: {PATHS_FILE}'))
source(PATHS_FILE, local = TRUE)
if (!exists('PATHS') || !is.list(PATHS)) stop(glue('config.R: {PATHS_FILE} must define a list named PATHS'))

`%||%` <- function(a, b) if (is.null(a)) b else a
expand_path <- function(p) if (is.null(p) || !nzchar(p)) NULL else path.expand(p)

input_path <- function(key) {
    p <- expand_path(PATHS[[key]])
    if (is.null(p)) stop(glue('scripts/paths.R: entry "{key}" is not set'))
    if (!file.exists(p)) stop(glue('scripts/paths.R: entry "{key}" points to a missing file: {p}'))
    p
}
optional_path <- function(key) { p <- expand_path(PATHS[[key]]); if (!is.null(p) && file.exists(p)) p else NULL }
optional_dir  <- function(key) { p <- optional_path(key); if (is.null(p)) NULL else sub('/*$', '/', p) }

OUTPUT_DIR <- sub('/*$', '/', expand_path(PATHS$output_dir) %||% stop('scripts/paths.R: output_dir is not set'))
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Analysis constants ----------------------------------------------------------------
SPLIT_SEED    <- 123
TRAIN_FRAC    <- 0.7
N_SUBTYPES    <- 6
KMEANS_NSTART <- 25
N_PGS_PCS     <- 8
N_BOOTSTRAP   <- 1000
FDR_STAGE1    <- 0.5
FDR_STAGE2    <- 0.1

ID_COL     <- 'X.projid'
COVARIATES <- c(paste0('PC', 1:10, '_AVG'), 'msex', 'age_death', 'age_bl', 'kronos')

# APOE gradient coding: ε2ε2 = -2 ... ε4ε4 = +2.
APOE_LEVELS <- c(`2` = 'ε4ε4', `1` = 'ε3ε4', `0` = 'ε3ε3', `-0.5` = 'ε2ε4', `-1` = 'ε2ε3', `-2` = 'ε2ε2')

STAGE1_PHENOTYPES <- c('cogn_global_lv', 'cogng_path_slope', 'amyloid', 'tangles', 'gpath')
STAGE2_PHENOTYPES <- c('plaq_n', 'plaq_d', 'plaq_n_ag', 'caa_4gp', 'gpath', 'cogdx', 'cognse_demog_slope', 'nft',
                       'tdp_st4', 'plaq_d_ag', 'nft_ag', 'cognwo_demog_slope', 'cogn_ps_lv', 'cognps_demog_slope',
                       'cogn_wo_lv', 'amyloid', 'cogn_se_lv', 'plaq_n_ec', 'nft_hip', 'nft_mt', 'plaq_n_hip', 'ceradsc',
                       'cogn_ep_lv', 'tangles', 'cognep_demog_slope', 'cogn_global_lv', 'gpath_3neocort',
                       'cogng_demog_slope', 'plaq_n_mf', 'plaq_d_mt', 'braaksc', 'cogn_po_lv', 'plaq_n_mt', 'nft_mf',
                       'nft_ec', 'plaq_d_mf')
PREDICTION_PHENOTYPES <- c('amyloid', 'nft', 'plaq_n', 'plaq_d', 'tangles', 'gpath',
                           'cogn_global_lv', 'niareagansc', 'cogdx')

# Figure 3: six phenotypes (x axis) and five PGS (rows), by display name; ids in the same order.
FIG3_PHENOTYPES <- c('Global cognitive function (19 tests)',
                     'Global AD pathology burden',
                     'Neuritic plaque burden (5 regions)',
                     'Amyloid level (% cortex area, 8 brain regions)',
                     'Diffuse plaque burden (5 regions)',
                     'Tangle density (IHC, 8 brain regions)')
FIG3_PGS     <- c("Alzheimer's/dementia (FH)", 'Apolipoprotein B', 'Use of sun/uv protection (Always)',
                  'Doctor diagnosed hayfever or allergic rhinitis', 'Prostate cancer')
FIG3_PGS_IDS <- c('FH1263', 'INI30640', 'BIN_FC40002267', 'BIN22126', 'cancer1044')

FIG2_SCATTER_PGS  <- c('INI30640_SCORE1_AVG', 'cancer1044_SCORE1_AVG')
FIG2_SCATTER_PHEN <- c('plaq_d', 'plaq_n')
GREAT_PGS <- c(INI30640 = 'Apolipoprotein B', cancer1044 = 'Prostate cancer')

# ---- House style (Figures 4 and 5): 7 pt black Helvetica, 0.75 pt strokes -------------------
HOUSE_FONT_FAMILY <- 'Helvetica'
HOUSE_FONT_PT     <- 7
HOUSE_LINE_PT     <- 0.75
HOUSE_LINE_MM     <- HOUSE_LINE_PT / ggplot2::.pt
HOUSE_TEXT_COLOUR <- 'black'
PANEL_IN <- list(
    a = c(2.05, 1.90), b = c(4.00, 1.90),
    c = c(2.15, 1.75), d = c(2.15, 1.75), e = c(2.15, 1.75),
    f = c(2.20, 2.10), g = c(2.60, 2.10), h = c(1.65, 2.10),
    prediction = c(7.10, 2.40))

theme_house <- function(base_size = HOUSE_FONT_PT, base_family = HOUSE_FONT_FAMILY) {
    ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
        text             = ggplot2::element_text(size = base_size, colour = HOUSE_TEXT_COLOUR, family = base_family),
        axis.text        = ggplot2::element_text(size = base_size, colour = HOUSE_TEXT_COLOUR),
        axis.title       = ggplot2::element_text(size = base_size, colour = HOUSE_TEXT_COLOUR),
        legend.text      = ggplot2::element_text(size = base_size, colour = HOUSE_TEXT_COLOUR),
        legend.title     = ggplot2::element_text(size = base_size, colour = HOUSE_TEXT_COLOUR),
        plot.title       = ggplot2::element_text(size = base_size, colour = HOUSE_TEXT_COLOUR),
        strip.text       = ggplot2::element_text(size = base_size, colour = HOUSE_TEXT_COLOUR),
        strip.background = ggplot2::element_blank(),
        axis.ticks       = ggplot2::element_line(linewidth = HOUSE_LINE_MM, colour = HOUSE_TEXT_COLOUR),
        axis.line        = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(linewidth = HOUSE_LINE_MM, colour = HOUSE_TEXT_COLOUR, fill = NA),
        panel.grid.major = ggplot2::element_line(linewidth = HOUSE_LINE_MM / 2, colour = 'grey85'),
        panel.grid.minor = ggplot2::element_blank(),
        legend.key.size  = grid::unit(base_size, 'pt'),
        legend.background = ggplot2::element_blank(),
        plot.margin      = grid::unit(c(2, 2, 2, 2), 'pt'))
}

# ---- Output helpers --------------------------------------------------------------------
open_pdf <- function(filename, width, height, ...) {
    if (capabilities('aqua')) grDevices::quartz(type = 'pdf', file = filename, width = width, height = height,
                                                family = HOUSE_FONT_FAMILY)
    else grDevices::cairo_pdf(filename, width = width, height = height, family = HOUSE_FONT_FAMILY)
}
open_panel_pdf <- function(path, panel) open_pdf(path, PANEL_IN[[panel]][1], PANEL_IN[[panel]][2])

fig_path <- function(figure, filename) {
    dir <- glue('{OUTPUT_DIR}{if (grepl("^FigureS", figure)) "supplementary" else "main"}/{figure}/')
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    glue('{dir}{filename}')
}
table_path <- function(filename) {
    dir <- glue('{OUTPUT_DIR}tables/')
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    glue('{dir}{filename}')
}

manifest_row <- function(figure, panel, path, status = 'generated', note = '') {
    data.frame(figure = figure, panel = panel, path = as.character(path), status = status,
               note = as.character(note), stringsAsFactors = FALSE)
}

say <- function(..., .envir = parent.frame()) cat(format(Sys.time(), '%H:%M:%S'), '|', glue(..., .envir = .envir), '\n')

options(device = 'pdf')
