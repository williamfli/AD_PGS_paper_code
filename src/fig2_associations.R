# src/fig2_associations.R -- the two-stage association and Figure 2.
#   make_associations()  Stage 1 (713 PGS x 5 global phenotypes), Stage 2 (12 PGS x 36 AD
#                        phenotypes), Stage 2 without the APOE region, APOE gradient alone;
#                        matrices under OUTPUT_DIR/association/, long format as Tables S3-S6
#   make_figure2()       Fig 2a-b heatmaps, Fig 2c-f scatter plots
#   make_figureS7()      Fig S7 Stage-1 heatmap of all 713 PGS
suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
})
ht_opt$message <- FALSE

# ======================= PART 1: the two-stage association ============================
# Every (predictor, phenotype) pair: mean-impute, z-score, OLS with covariates, BH FDR over all pairs.
fit_association_matrices <- function(data, predictors, phenotypes, covariates = COVARIATES) {
    stopifnot(all(c(predictors, phenotypes, covariates) %in% colnames(data)))
    for (column in c(predictors, phenotypes, covariates)) {
        v <- data[[column]]
        if (anyNA(v)) v[is.na(v)] <- mean(v, na.rm = TRUE)
        data[[column]] <- v
    }
    zscore <- function(v) (v - mean(v)) / sd(v)
    X_cov <- as.matrix(data[, covariates])
    n <- nrow(data)
    p <- ncol(X_cov) + 2
    Z_pred <- vapply(predictors, function(v) zscore(data[[v]]), numeric(n))
    Z_phen <- vapply(phenotypes, function(v) zscore(data[[v]]), numeric(n))

    effects <- errors <- pvals <- matrix(NA_real_, length(predictors), length(phenotypes),
                                         dimnames = list(predictors, phenotypes))
    for (i in seq_along(predictors)) {
        X <- cbind(`(Intercept)` = 1, x = Z_pred[, i], X_cov)
        qrX <- qr(X)
        XtX_inv_xx <- chol2inv(qr.R(qrX))[2, 2]
        coef <- qr.coef(qrX, Z_phen)
        resid <- Z_phen - X %*% coef
        sigma2 <- colSums(resid^2) / (n - p)
        se <- sqrt(XtX_inv_xx * sigma2)
        effects[i, ] <- coef[2, ]
        errors[i, ] <- se
        pvals[i, ] <- 2 * pt(abs(coef[2, ] / se), df = n - p, lower.tail = FALSE)
    }
    fdr <- matrix(p.adjust(as.vector(t(pvals)), method = 'fdr'), nrow(pvals), ncol(pvals),
                  byrow = TRUE, dimnames = dimnames(pvals))
    list(effects = as.data.frame(effects), errors = as.data.frame(errors),
         pvals = as.data.frame(pvals), fdr = as.data.frame(fdr), n = n)
}

read_scores <- function(key) {
    d <- read.delim(input_path(key), sep = '\t', header = TRUE, check.names = FALSE, colClasses = c('#projid' = 'character'))
    d[['#projid']] <- formatC(as.integer(d[['#projid']]), width = 8, flag = '0', format = 'd')
    d
}

cohort_with_scores <- function(scores, suffix) {
    cohort <- load_cohort()
    score_cols <- setdiff(colnames(scores), '#projid')
    idx <- match(cohort[[ID_COL]], scores[['#projid']])
    if (anyNA(idx)) stop(glue('{sum(is.na(idx))} cohort individuals missing from the score table'))
    add <- scores[idx, score_cols, drop = FALSE]
    colnames(add) <- paste0(score_cols, suffix)
    list(data = cbind(cohort, add), predictors = colnames(add))
}

write_association_run <- function(result, run, predictor_names = NULL) {
    dir <- glue('{OUTPUT_DIR}association/{run}/')
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    for (name in c('effects', 'errors', 'pvals', 'fdr')) {
        m <- result[[name]]
        if (!is.null(predictor_names)) rownames(m) <- predictor_names
        write.table(m, glue('{dir}{name}.tsv'), sep = '\t')
    }
    dir
}

# One row per (predictor, phenotype) pair with display names, for Tables S3-S6.
association_long_table <- function(result, predictor_names = NULL, feature_column = 'UKB PGS',
                                   feature_labels = mapper_pgs, sort_by_feature = FALSE) {
    cell <- function(m) { m <- as.matrix(m); if (!is.null(predictor_names)) rownames(m) <- predictor_names; m }
    e <- cell(result$effects)
    g <- expand.grid(phenotype = colnames(e), predictor = rownames(e), stringsAsFactors = FALSE)
    pick <- function(m) cell(m)[cbind(g$predictor, g$phenotype)]
    out <- data.frame(feature = unname(feature_labels(g$predictor)), phenotype = unname(mapper_pheno(g$phenotype)),
                      `Effect Size` = pick(result$effects), `Standard Error` = pick(result$errors),
                      `P-Value` = pick(result$pvals), FDR = pick(result$fdr), check.names = FALSE, stringsAsFactors = FALSE)
    names(out)[1:2] <- c(feature_column, 'Observed ROSMAP Phenotype')
    if (sort_by_feature) out <- out[order(tolower(out[[1]]), out$FDR, method = 'radix'), ]
    rownames(out) <- NULL
    out
}

make_associations <- function() {
    cohort <- load_cohort()
    pgs_cols <- grep('_SCORE1_AVG$', colnames(cohort), value = TRUE)
    stopifnot(length(pgs_cols) == 12)
    rows <- list()

    run <- function(key, data, predictors, phenotypes, predictor_names = NULL, label) {
        say('Association {key}: {length(predictors)} predictors x {length(phenotypes)} phenotypes')
        res <- fit_association_matrices(data, predictors, phenotypes)
        dir <- write_association_run(res, key, predictor_names)
        rows[[key]] <<- manifest_row('Associations', key, glue('{dir}{{effects,errors,pvals,fdr}}.tsv'), note = glue('{label}; n = {res$n}'))
        res
    }

    s2   <- run('stage2', cohort, pgs_cols, STAGE2_PHENOTYPES, label = '12 PGS x 36 phenotypes')
    apoe <- run('apoe_only', cohort, 'apoe_genotype', STAGE2_PHENOTYPES, label = 'APOE gradient x 36 phenotypes')

    na <- cohort_with_scores(read_scores('pgs_scores_no_apoe'), suffix = '__noAPOE')
    na_names <- sub('__noAPOE$', '', na$predictors)
    s2na <- run('stage2_no_apoe', na$data, na$predictors, STAGE2_PHENOTYPES, na_names,
                label = '12 APOE-excluded PGS x 36 phenotypes')

    s1in <- cohort_with_scores(read_scores('all_pgs_scores'), suffix = '__all')
    s1_names <- sub('__all$', '', s1in$predictors)
    s1 <- run('stage1', s1in$data, s1in$predictors, STAGE1_PHENOTYPES, s1_names, label = '713 PGS x 5 global phenotypes')

    # The Stage-1 rule (any FDR < 0.5) must select exactly the 12 PGS carried in the cohort table.
    prioritized <- s1_names[apply(s1$fdr, 1, function(r) any(r < FDR_STAGE1))]
    if (!setequal(prioritized, pgs_cols))
        stop(glue('Stage-1 prioritization selects {length(prioritized)} PGS, not the 12 in the cohort table: ',
                  '{paste(setdiff(prioritized, pgs_cols), collapse = ", ")}'))
    say('  Stage-1 rule reproduces the 12 prioritized PGS; Stage-2 associations at FDR < {FDR_STAGE2}: {sum(s2$fdr < FDR_STAGE2)}')

    tables <- list(
        list('TableS3', 'TableS3_stage1_associations.tsv',         association_long_table(s1, s1_names, sort_by_feature = TRUE)),
        list('TableS4', 'TableS4_stage2_associations.tsv',         association_long_table(s2)),
        list('TableS5', 'TableS5_stage2_associations_no_apoe.tsv', association_long_table(s2na, na_names)),
        list('TableS6', 'TableS6_apoe_gradient_associations.tsv',  association_long_table(apoe, feature_column = 'Genetic Feature',
                                                                                          feature_labels = function(x) rep('APOE-gradient', length(x)))))
    for (t in tables) {
        path <- table_path(t[[2]])
        write.table(t[[3]], path, sep = '\t', quote = FALSE, row.names = FALSE)
        rows[[t[[1]]]] <- manifest_row(t[[1]], 'table', path, note = glue('{nrow(t[[3]])} rows'))
    }
    do.call(rbind, rows)
}

# ======================= PART 2: Figure 2 and Figure S7 ===============================
ASSOC_COL_FUN <- function() colorRamp2(c(-0.4, 0, 0.4), c('blue', 'white', 'orange'))

significance_marks <- function(fdr, thresholds = c(0.05, 0.01, 0.001)) {
    fdr <- as.matrix(fdr)
    marks <- matrix('', nrow = nrow(fdr), ncol = ncol(fdr), dimnames = dimnames(fdr))
    for (k in seq_along(thresholds)) marks[which(fdr < thresholds[k])] <- strrep('*', k)
    marks
}

significance_legend <- function(thresholds = c(0.05, 0.01, 0.001), fontsize = 9, title = 'Significance of association') {
    k <- rev(seq_along(thresholds))
    Legend(title = title,
           labels = sprintf('FDR < %s', format(thresholds[k], scientific = FALSE, drop0trailing = TRUE)),
           graphics = lapply(k, function(i) {
               force(i)
               function(x, y, w, h) grid.text(strrep('*', i), x, y, gp = gpar(fontsize = fontsize))
           }),
           title_gp = gpar(fontsize = fontsize, fontface = 'plain'), labels_gp = gpar(fontsize = fontsize))
}

marks_cell_fun <- function(marks, fontsize) {
    function(j, i, x, y, width, height, fill) {
        if (marks[i, j] != '') grid.text(marks[i, j], x, y, gp = gpar(fontsize = fontsize))
    }
}

# Title of each slice of a split heatmap: the name in `anchors` whose dimname falls in that slice.
slice_titles <- function(orders, dimnames, anchors, default = '') {
    vapply(orders, function(idx) {
        hit <- names(anchors)[anchors %in% dimnames[idx]]
        if (length(hit) == 1) hit else default
    }, character(1))
}

drawn_size <- function(draw_fun, pad_in = 0.3) {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    ht <- draw_fun()
    c(width = convertWidth(ComplexHeatmap:::width(ht), 'inch', valueOnly = TRUE) + 2 * pad_in,
      height = convertHeight(ComplexHeatmap:::height(ht), 'inch', valueOnly = TRUE) + 2 * pad_in)
}

# ---- Figure 2a/b ------------------------------------------------------------------------
# Row (phenotype) and column (PGS) dendrograms of the Stage-2 matrix; also used by Figure S8.
stage2_dendrograms <- function() {
    m <- t(as.matrix(load_stage2()$effects))
    set.seed(123)
    probe <- prepare(Heatmap(m, cluster_rows = TRUE, cluster_columns = TRUE))
    suppressWarnings(list(row = row_dend(probe), col = column_dend(probe)))
}
stage2_pgs_dendrogram <- function() rev(stage2_dendrograms()$col)

build_figure2_heatmaps <- function(cell_mm = 5.5, fontsize = 9) {
    s1 <- load_stage1_prioritized()
    s2 <- load_stage2()
    m1 <- t(as.matrix(s1$effects)); k1 <- significance_marks(t(s1$fdr))
    m2 <- t(as.matrix(s2$effects)); k2 <- significance_marks(t(s2$fdr))
    stopifnot(setequal(colnames(m1), colnames(m2)))
    col_dend <- stage2_pgs_dendrogram()
    m1 <- m1[, colnames(m2)]; k1 <- k1[rownames(m1), colnames(m1)]
    k2 <- k2[rownames(m2), colnames(m2)]

    cell <- unit(cell_mm, 'mm')
    legend_param <- list(title = 'Strength of UKB PGS association with\nROSMAP AD phenotype',
                         legend_direction = 'horizontal', legend_width = unit(6.5, 'cm'),
                         at = c(-0.4, -0.2, 0, 0.2, 0.4),
                         title_gp = gpar(fontsize = fontsize, fontface = 'plain'),
                         labels_gp = gpar(fontsize = fontsize), title_position = 'topcenter')

    col_anchors <- c('AD discordant' = 'Prostate cancer', 'AD concordant' = 'Apolipoprotein B', 'AD' = "Alzheimer's/dementia (FH)")
    set.seed(123)
    probe <- prepare(Heatmap(m2, cluster_rows = TRUE, row_split = 3, cluster_columns = col_dend, column_split = 3))
    col_titles <- slice_titles(column_order(probe), colnames(m2), col_anchors)

    common_ht <- list(col = ASSOC_COL_FUN(), border = TRUE,
                      row_names_side = 'left', row_names_gp = gpar(fontsize = fontsize),
                      row_names_max_width = unit(12, 'cm'), column_names_max_height = unit(12, 'cm'),
                      row_dend_side = 'right', row_dend_width = unit(8, 'mm'),
                      row_title_side = 'left', row_title_gp = gpar(fontsize = fontsize + 1),
                      heatmap_legend_param = legend_param)
    shared_columns <- list(cluster_columns = col_dend, column_split = 3, column_gap = unit(0, 'mm'),
                           column_dend_height = unit(12, 'mm'), cluster_column_slices = TRUE,
                           column_title = col_titles, column_title_gp = gpar(fontsize = fontsize + 1),
                           column_names_gp = gpar(fontsize = fontsize))

    ht_a <- do.call(Heatmap, c(
        list(m1, name = 'effect', width = ncol(m1) * cell, height = nrow(m1) * cell, cluster_rows = TRUE,
             row_title = 'Stage 1', cell_fun = marks_cell_fun(k1, fontsize - 2)),
        shared_columns, common_ht))
    ht_b <- do.call(Heatmap, c(
        list(m2, name = 'effect_stage2', width = ncol(m2) * cell, height = nrow(m2) * cell,
             cluster_rows = TRUE, row_split = 3, row_gap = unit(0, 'mm'), row_title = 'Stage 2',
             cell_fun = marks_cell_fun(k2, fontsize - 2)),
        shared_columns, common_ht))

    list(ht_a = ht_a, ht_b = ht_b, fontsize = fontsize, sig_legend = significance_legend(fontsize = fontsize))
}

draw_figure2_heatmaps <- function(built, which = c('a', 'b')) {
    which <- match.arg(which)
    fs <- built$fontsize
    draw(if (which == 'a') built$ht_a else built$ht_b, main_heatmap = 1,
         heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom', merge_legend = TRUE,
         annotation_legend_list = list(built$sig_legend), legend_gap = unit(12, 'mm'),
         row_title = 'Observed ROSMAP AD phenotypes', row_title_side = 'left', row_title_gp = gpar(fontsize = fs + 1),
         column_title = '12 prioritized AD-relevant UKB PGS', column_title_side = 'bottom', column_title_gp = gpar(fontsize = fs + 1),
         padding = unit(c(2, 2, 2, 2), 'mm'))
}

# ---- Figure 2c-f: scatter plots (line = Stage-2 beta x SD(phenotype), ribbon = 1.96 SE) -------
draw_figure2_scatter <- function(merged_data, pgs_vars, phenotype_vars, effects, errors, fdr) {
    stopifnot(length(pgs_vars) == 2, length(phenotype_vars) == 2)
    op <- par(mfrow = c(2, 2), mar = c(4.2, 4.5, 1.5, 1), oma = c(2.5, 0, 0, 0), mgp = c(2.6, 0.8, 0))
    on.exit(par(op), add = TRUE)

    standardized <- function(pgs, pheno) {
        x <- merged_data[[pgs]]; y <- merged_data[[pheno]]
        ok <- !is.na(x) & !is.na(y)
        x <- x[ok]; y <- y[ok]
        list(x = (x - mean(x)) / sd(x), y = y)
    }
    x_lim <- lapply(pgs_vars, function(p) range(unlist(lapply(phenotype_vars, function(q) standardized(p, q)$x))))
    y_lim <- lapply(phenotype_vars, function(q) range(unlist(lapply(pgs_vars, function(p) standardized(p, q)$y))))

    for (j in seq_along(phenotype_vars)) {
        for (i in seq_along(pgs_vars)) {
            d <- standardized(pgs_vars[i], phenotype_vars[j])
            pgs_label <- mapper_pgs(pgs_vars[i]); pheno_label <- mapper_pheno(phenotype_vars[j])
            beta <- as.numeric(effects[pgs_label, pheno_label])
            se   <- as.numeric(errors[pgs_label, pheno_label])
            q    <- as.numeric(fdr[pgs_label, pheno_label])
            plot(d$x, d$y, xlim = x_lim[[i]], ylim = y_lim[[j]], xlab = pgs_label, ylab = pheno_label,
                 pch = 19, col = rgb(0, 0, 1, 0.3), cex = 0.8, cex.lab = 1.05, cex.axis = 0.95, las = 1)
            xr <- x_lim[[i]]; sy <- sd(d$y); my <- mean(d$y)
            polygon(c(xr, rev(xr)), c((beta + 1.96 * se) * sy * xr + my, rev((beta - 1.96 * se) * sy * xr + my)),
                    col = rgb(1, 0, 0, 0.2), border = NA)
            lines(xr, beta * sy * xr + my, col = 'red', lwd = 3)
            legend('topleft', inset = c(0.02, 0.06), bty = 'n', cex = 0.95,
                   legend = c('Standardized regression', sprintf('coefficient: %.3f ± %.3f', beta, se), '',
                              'False discovery rate (q-value):', format(q, scientific = TRUE, digits = 2)))
        }
    }
    mtext('Standardized UK Biobank PGS (z-score)', side = 1, outer = TRUE, line = 1, cex = 1.05)
    invisible(NULL)
}

make_figure2 <- function() {
    say('Figure 2: association heatmaps')
    built <- build_figure2_heatmaps()
    rows <- list()

    for (spec in list(list(which = 'a', file = 'Fig2a_stage1_heatmap.pdf', panel = 'a', note = 'Stage-1 associations of the 12 prioritized PGS; shared Stage-2 PGS dendrogram'),
                      list(which = 'b', file = 'Fig2b_stage2_heatmap.pdf', panel = 'b', note = 'Stage-2 associations, 12 PGS x 36 phenotypes'))) {
        size <- drawn_size(function() draw_figure2_heatmaps(built, spec$which))
        path <- fig_path('Figure2', spec$file)
        pdf(path, width = size[['width']], height = size[['height']], useDingbats = FALSE)
        draw_figure2_heatmaps(built, spec$which)
        dev.off()
        rows[[length(rows) + 1]] <- manifest_row('Figure2', spec$panel, path, note = spec$note)
    }

    say('Figure 2: scatter plots')
    s2 <- load_stage2()
    path <- fig_path('Figure2', 'Fig2cdef_scatter_plots.pdf')
    pdf(path, width = 9, height = 9, useDingbats = FALSE)
    draw_figure2_scatter(load_cohort(), FIG2_SCATTER_PGS, FIG2_SCATTER_PHEN, s2$effects, s2$errors, s2$fdr)
    dev.off()
    rows[[length(rows) + 1]] <- manifest_row('Figure2', 'c-f', path, note = 'ApoB and prostate-cancer PGS vs diffuse/neuritic plaque burden')
    do.call(rbind, rows)
}

# ---- Figure S7: Stage-1 heatmap of all 713 PGS; '*' = FDR < 0.5 ------------------------------
make_figureS7 <- function() {
    say('Figure S7: Stage-1 heatmap of all PGS')
    s1 <- load_stage1()
    m <- t(as.matrix(s1$effects))
    marks <- significance_marks(t(s1$fdr), thresholds = FDR_STAGE1)
    fontsize <- 9
    set.seed(123)
    ht <- Heatmap(m, name = 'effect', col = ASSOC_COL_FUN(),
                  width = ncol(m) * unit(0.5, 'mm'), height = nrow(m) * unit(7, 'mm'),
                  cluster_rows = TRUE, cluster_columns = TRUE,
                  row_dend_side = 'left', row_dend_width = unit(6, 'mm'), column_dend_height = unit(15, 'mm'),
                  show_column_names = FALSE,
                  row_names_side = 'right', row_names_gp = gpar(fontsize = fontsize), row_names_max_width = unit(12, 'cm'),
                  row_title = 'Selected ROSMAP phenotypes', row_title_gp = gpar(fontsize = fontsize + 1),
                  column_title = glue('UK Biobank PGS (n = {ncol(m)})'), column_title_side = 'bottom',
                  column_title_gp = gpar(fontsize = fontsize + 1),
                  cell_fun = marks_cell_fun(marks, 3),
                  heatmap_legend_param = list(title = 'Effect size (beta)', legend_direction = 'horizontal',
                                              legend_width = unit(4, 'cm'), at = c(-0.4, -0.2, 0, 0.2, 0.4),
                                              title_gp = gpar(fontsize = fontsize, fontface = 'plain'),
                                              labels_gp = gpar(fontsize = fontsize)))
    sig_legend <- significance_legend(thresholds = FDR_STAGE1, fontsize = fontsize, title = 'Significance')
    draw_it <- function() draw(ht, heatmap_legend_side = 'bottom', annotation_legend_side = 'bottom', merge_legend = TRUE,
                               annotation_legend_list = list(sig_legend), legend_gap = unit(12, 'mm'),
                               padding = unit(c(2, 2, 2, 2), 'mm'))
    size <- drawn_size(draw_it)
    path <- fig_path('FigureS7', 'FigS7_stage1_all_713_pgs_heatmap.pdf')
    pdf(path, width = size[['width']], height = size[['height']], useDingbats = FALSE)
    draw_it()
    dev.off()
    manifest_row('FigureS7', 'S7', path, note = glue('Stage-1 associations of all {ncol(m)} PGS; * = FDR < {FDR_STAGE1}'))
}
