# src/fig5_pgs_subtypes.R -- Figure 5 and everything built on the PGS PCA and polygenic subtypes.
#   make_figure5_pca()             Fig 5a-e: PGS PCA fitted on the seeded 70% training split
#   make_figureS5()                Fig S5a-b: individual-level PGS and PGS-PC heatmaps
#   make_figure5_subtypes()        Fig 5f-h: subtype biplot, plaque violins (subtype 1 vs 6), APOE allelotypes
#   make_figureS11()               Fig S11a: subtype means of six standardized phenotypes; S11b-d: pairwise subtype tests
#   make_figureS12()               Fig S12: Fig 5g within ε3/ε3 individuals
#   make_figureS13_and_tables()    Fig S13 bootstrap examples, Table S11 stability
#   make_figureS14()               Fig S14: Fig 5g on the test split, k-means fit on the training split only

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(patchwork)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
    library(dendextend)
})

# ============================ PART 1: PGS PCA (Fig 5a-e, S5) ==========================
pc_axis_label <- function(k, variance_explained) sprintf('PGS PC%d (%.1f%%)', k, 100 * variance_explained[k])

pgs_loadings <- function(fit = fit_pgs_pca()) {
    rot <- fit$pca$rotation
    rownames(rot) <- sapply(rownames(rot), mapper_pgs)
    rot
}

# House style in grid units: lwd is 1/96 in, hence HOUSE_LINE_PT / 0.75.
fig5_gp <- function(...) gpar(fontsize = HOUSE_FONT_PT, col = HOUSE_TEXT_COLOUR, fontfamily = HOUSE_FONT_FAMILY, ...)
FIG5_LWD <- HOUSE_LINE_PT / 0.75
FIG5_TEXT_SIZE <- HOUSE_FONT_PT / ggplot2::.pt
FIG5_ARROW <- arrow(length = unit(1.5, 'mm'))

FIG5_ARROW_LABELS <- c(
    'Types of physical activity in last 4 weeks (other exercises)' = 'Physical activity\n(other exercises)',
    'Doctor diagnosed hayfever or allergic rhinitis'                = 'Hayfever or\nallergic rhinitis',
    'Vol. of grey matter in Planum Polare (R)'                      = 'Grey matter vol.\nPlanum Polare (R)',
    'Use of sun/uv protection (Always)'                             = 'Sun/UV protection\n(Always)',
    "AD alzheimer's disease"                                        = "AD alzheimer's\ndisease",
    "Alzheimer's/dementia (FH)"                                     = "Alzheimer's/\ndementia (FH)")
FIG5_HEATMAP_LABELS <- c('Types of physical activity in last 4 weeks (other exercises)' = 'Physical activity (other exercises)')
relabel <- function(x, map) ifelse(x %in% names(map), map[x], x)

# Arrows shorter than this fraction of the panel's longest arrow are drawn but not labelled.
FIG5_ARROW_LABEL_MIN_FRAC <- c(c = 0.25, d = 0.25, e = 0.33)

# ---- Figure 5a: variance explained --------------------------------------------------------
plot_variance_explained <- function(fit) {
    ve <- 100 * fit$variance_explained
    df <- data.frame(PC = factor(seq_along(ve), levels = seq_along(ve)), variance = as.numeric(ve), cumulative = cumsum(as.numeric(ve)))
    df$bar_label <- ifelse(df$variance == df$cumulative, '', sprintf('%.1f', df$variance))
    ggplot(df, aes(x = PC, y = variance)) +
        geom_col(fill = '#d9642a', color = '#B22222', linewidth = HOUSE_LINE_MM, width = 0.85) +
        geom_text(aes(y = variance + 2, label = bar_label), angle = 90, hjust = 0, vjust = 0.5,
                  size = FIG5_TEXT_SIZE, family = HOUSE_FONT_FAMILY, color = '#d9642a') +
        geom_line(aes(y = cumulative, group = 1), color = '#B22222', linewidth = HOUSE_LINE_MM) +
        geom_point(aes(y = cumulative), color = '#B22222', size = 0.8) +
        geom_text(aes(y = cumulative + 3, label = sprintf('%.1f', cumulative)), angle = 90, hjust = 0, vjust = 0.5,
                  size = FIG5_TEXT_SIZE, family = HOUSE_FONT_FAMILY, color = '#B22222') +
        annotate('text', x = 12.5, y = 62, hjust = 1, label = 'Cumulative\nvariance (%)', lineheight = 0.9,
                 color = '#B22222', size = FIG5_TEXT_SIZE, family = HOUSE_FONT_FAMILY) +
        annotate('text', x = 12.5, y = 36, hjust = 1, label = 'Individual\ncomponents (%)', lineheight = 0.9,
                 color = '#d9642a', size = FIG5_TEXT_SIZE, family = HOUSE_FONT_FAMILY) +
        scale_y_continuous(limits = c(0, 128), breaks = seq(0, 100, 25), expand = c(0, 0)) +
        scale_x_discrete(expand = expansion(add = 0.6)) +
        labs(x = 'PGS principal components (PCs)', y = 'Proportion of variance explained (%)') +
        theme_house()
}

# ---- Figure 5b: |loading| heatmap with the published dendrogram leaf order ---------------------
loadings_heatmap_dendrograms <- function(mat) {
    dd_cols <- as.dendrogram(hclust(dist(t(mat)), method = 'ward.D2'))
    dd_cols <- reorder(dd_cols, c(1, 4, 5, 9, 11, 10, 12, 7, 6, 3, 8, 2))
    dd_cols <- flip_leaves(dd_cols, c(10), c(1, 12))
    dd_cols <- flip_leaves(dd_cols, c(9), c(2, 3))
    dd_rows <- reorder(as.dendrogram(hclust(dist(mat), method = 'ward.D2')), c(4, 6, 5, 7, 10, 9, 12, 3, 2, 8, 1, 11))
    list(rows = dd_rows, cols = dd_cols)
}

draw_loadings_heatmap <- function(fit) {
    mat <- abs(pgs_loadings(fit))
    dd <- loadings_heatmap_dendrograms(mat)
    rownames(mat) <- relabel(rownames(mat), FIG5_HEATMAP_LABELS)
    col_fun <- colorRamp2(c(0, max(mat)), c('white', 'blue'))
    ht <- Heatmap(mat, name = 'Loading', col = col_fun, show_heatmap_legend = FALSE,
                  row_names_side = 'left', column_names_side = 'bottom', column_names_rot = 90,
                  column_title = 'PGS principal components (PCs)', column_title_side = 'bottom',
                  row_title = 'Prioritized UKB PGS',
                  row_names_gp = fig5_gp(), column_names_gp = fig5_gp(), row_title_gp = fig5_gp(), column_title_gp = fig5_gp(),
                  row_dend_side = 'right', cluster_rows = dd$rows, cluster_columns = dd$cols,
                  row_dend_width = unit(0.15, 'in'), column_dend_height = unit(0.32, 'in'),
                  row_dend_gp = gpar(lwd = FIG5_LWD, col = HOUSE_TEXT_COLOUR),
                  column_dend_gp = gpar(lwd = FIG5_LWD, col = HOUSE_TEXT_COLOUR),
                  border = TRUE, border_gp = gpar(lwd = FIG5_LWD, col = HOUSE_TEXT_COLOUR),
                  row_names_max_width = unit(2.2, 'in'), column_names_max_height = unit(0.3, 'in'))
    legend <- Legend(col_fun = col_fun, at = c(0, 0.5, 1), labels = c('0', '0.5', '1'),
                     title = 'PGS loadings on principal components', direction = 'horizontal',
                     title_gp = fig5_gp(), labels_gp = fig5_gp(), title_position = 'topleft',
                     legend_width = unit(0.6, 'in'), grid_height = unit(2, 'mm'))
    draw(ht, padding = unit(c(1, 1, 1, 1), 'mm'))
    draw(legend, x = unit(1, 'mm'), y = unit(1, 'npc') - unit(1, 'mm'), just = c('left', 'top'))
}

# ---- Figure 5c-e: test-split biplots coloured by neuritic plaque burden (NA dropped) ---------
plot_pgs_biplot <- function(fit, pc_x, pc_y, arrow_scale, colour_pheno = 'plaq_n', legend = TRUE, label_min_frac = 0.25) {
    cohort <- load_cohort()
    scores <- as.data.frame(fit$scores_all[fit$test_idx, , drop = FALSE])
    scores$pheno <- cohort[fit$test_idx, colour_pheno]
    scores <- scores[!is.na(scores$pheno), ]
    loadings <- as.data.frame(pgs_loadings(fit))
    loadings$var <- relabel(rownames(loadings), FIG5_ARROW_LABELS)
    x <- paste0('PC', pc_x); y <- paste0('PC', pc_y)
    loadings$xend <- arrow_scale * loadings[[x]]; loadings$yend <- arrow_scale * loadings[[y]]
    arrow_len <- sqrt(loadings$xend^2 + loadings$yend^2)
    labelled <- loadings[arrow_len >= label_min_frac * max(arrow_len), ]
    range_x <- max(abs(c(scores[[x]], loadings$xend)))
    range_y <- max(abs(c(scores[[y]], loadings$yend)))
    p <- ggplot(scores, aes(x = .data[[x]], y = .data[[y]], color = pheno)) +
        geom_point(size = 0.8) +
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = xend, yend = yend),
                     arrow = FIG5_ARROW, linewidth = HOUSE_LINE_MM, color = 'blue', inherit.aes = FALSE) +
        geom_text_repel(data = labelled, aes(x = xend, y = yend, label = var),
                        size = FIG5_TEXT_SIZE, family = HOUSE_FONT_FAMILY, color = HOUSE_TEXT_COLOUR,
                        lineheight = 0.85, segment.size = HOUSE_LINE_MM, segment.color = HOUSE_TEXT_COLOUR,
                        box.padding = 0.2, point.padding = 0.1, min.segment.length = 0,
                        force = 5, force_pull = 0.2, max.iter = 1e6, max.time = 30,
                        max.overlaps = Inf, seed = SPLIT_SEED, inherit.aes = FALSE) +
        scale_color_gradient(low = 'white', high = 'red', name = 'Neuritic plaque burden') +
        labs(x = pc_axis_label(pc_x, fit$variance_explained), y = pc_axis_label(pc_y, fit$variance_explained)) +
        theme_house() +
        coord_cartesian(xlim = c(-range_x, range_x), ylim = c(-range_y, range_y), clip = 'off')
    if (legend) {
        p + theme(legend.position = 'top', legend.direction = 'horizontal', legend.title.position = 'left',
                  legend.key.height = unit(2, 'mm'), legend.key.width = unit(2.4, 'mm'),
                  legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, 'mm'),
                  legend.ticks = element_line(linewidth = HOUSE_LINE_MM, colour = HOUSE_TEXT_COLOUR),
                  legend.frame = element_rect(linewidth = HOUSE_LINE_MM, colour = HOUSE_TEXT_COLOUR))
    } else {
        p + theme(legend.position = 'none')
    }
}

make_figure5_pca <- function() {
    say('Figure 5a-e: PGS PCA on the training split')
    fit <- fit_pgs_pca()
    rows <- list()

    path <- fig_path('Figure5', 'Fig5a_variance_explained.pdf')
    open_panel_pdf(path, 'a'); print(plot_variance_explained(fit)); dev.off()
    rows[[1]] <- manifest_row('Figure5', 'a', path, note = 'variance explained by the 12 PGS PCs (PCA on the 70% training split)')

    path <- fig_path('Figure5', 'Fig5b_pgs_loadings_heatmap.pdf')
    open_panel_pdf(path, 'b'); draw_loadings_heatmap(fit); dev.off()
    rows[[2]] <- manifest_row('Figure5', 'b', path, note = 'absolute PCA loadings, 12 PGS x 12 PCs, published dendrogram leaf order')

    panels <- list(c = c(1, 2, 7), d = c(1, 6, 5), e = c(2, 6, 5))
    n_plotted <- sum(!is.na(load_cohort()[fit$test_idx, 'plaq_n']))
    for (panel in names(panels)) {
        p <- panels[[panel]]
        path <- fig_path('Figure5', sprintf('Fig5%s_biplot_PC%d_PC%d.pdf', panel, p[1], p[2]))
        open_panel_pdf(path, panel)
        print(plot_pgs_biplot(fit, p[1], p[2], arrow_scale = p[3], legend = (panel == 'c'), label_min_frac = FIG5_ARROW_LABEL_MIN_FRAC[[panel]]))
        dev.off()
        rows[[length(rows) + 1]] <- manifest_row('Figure5', panel, path, note = sprintf(
            'test-split individuals with plaq_n (n = %d of %d) on PC%d vs PC%d, coloured by neuritic plaque burden; loadings x%d',
            n_plotted, length(fit$test_idx), p[1], p[2], p[3]))
    }
    do.call(rbind, rows)
}

# ---- Figure S5: individual-level heatmaps ----------------------------------------------------
z_columns <- function(m) apply(m, 2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

individuals_heatmap <- function(mat, col_fun, row_title, legend_title) {
    Heatmap(mat, width = ncol(mat) * unit(0.24, 'mm'), height = nrow(mat) * unit(7, 'mm'),
            column_dend_height = unit(2, 'cm'), row_dend_width = unit(2, 'cm'),
            col = col_fun, row_title = row_title, column_title = 'Individuals',
            row_title_gp = gpar(fontsize = 22), column_title_gp = gpar(fontsize = 22),
            cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = TRUE, show_column_names = FALSE,
            row_names_gp = gpar(fontsize = 22),
            heatmap_legend_param = list(title = legend_title, legend_direction = 'horizontal',
                                        title_gp = gpar(fontsize = 22), labels_gp = gpar(fontsize = 22)))
}

make_figureS5 <- function() {
    say('Figure S5a-b: PGS and PGS-PC heatmaps across individuals')
    fit <- fit_pgs_pca()
    cohort <- load_cohort()

    pgs_mat <- t(z_columns(as.matrix(cohort[, fit$pgs_cols])))
    rownames(pgs_mat) <- sapply(rownames(pgs_mat), mapper_pgs)
    colnames(pgs_mat) <- cohort[[ID_COL]]
    path_a <- fig_path('FigureS5', 'FigS5a_pgs_heatmap.pdf')
    pdf(path_a, width = 35, height = 25)
    draw(individuals_heatmap(pgs_mat, colorRamp2(c(-2, 0, 2), c('navy', 'white', 'gold')), 'Polygenic scores', 'Z-score'),
         heatmap_legend_side = 'top', padding = unit(c(4, 4, 4, 4), 'cm'))
    dev.off()

    pc_mat <- t(fit$scores_all)
    colnames(pc_mat) <- cohort[[ID_COL]]
    path_b <- fig_path('FigureS5', 'FigS5b_pgs_pc_heatmap.pdf')
    pdf(path_b, width = 35, height = 25)
    draw(individuals_heatmap(pc_mat, colorRamp2(c(-2, 0, 2), c('#0072B2', 'white', '#D55E00')),
                             'PGS principal components (PCs)', 'Principal component (PC)'),
         heatmap_legend_side = 'top', padding = unit(c(4, 4, 4, 4), 'cm'))
    dev.off()

    rbind(manifest_row('FigureS5', 'a', path_a, note = glue('12 prioritized PGS z-scored over {ncol(pgs_mat)} individuals, both axes clustered')),
          manifest_row('FigureS5', 'b', path_b, note = 'PGS PC1-12 scores of all individuals, both axes clustered'))
}

# ============== PART 2: subtypes (Fig 5f-h, S11-S13, Table S11) =======================
SUBTYPE_PHENOTYPE_IDS <- c('cogn_global_lv', 'gpath', 'plaq_n', 'amyloid', 'plaq_d', 'tangles')
PLAQUE_IDS <- c(plaq_d = 'Diffuse plaque', plaq_n = 'Neuritic plaque')
VIOLIN_SUBTYPES <- c(1, 6)   # the two subtypes compared in Fig 5g, S11 (boxed cells), S12 and S14
PAIRWISE_PHENOTYPE_IDS <- c('plaq_n', 'plaq_d', 'nft')   # Fig S11b-d
PAIRWISE_DIFF_CAP <- 0.8   # colour limit (+/-) for differences in subtype mean z-score
PAIRWISE_BAR_FRACTION <- 0.158   # colour bar inset height as a fraction of the Fig S11b-d page
S11_PAIRWISE_IN <- c(4.15, 5.7)    # Fig S11b-d page (in)
LOADING_ARROW_SCALE <- 7
S13_EXAMPLE_ITERATIONS <- c(1, 27, 77)
UNLABELLED_PGS_PATTERN <- 'physical activity|sun/uv|grey matter'

subtype_colours <- function() setNames(scales::hue_pal()(6)[c(5, 1, 6, 2, 4, 3)], 1:6)

# ---- The subtype pipeline ---------------------------------------------------------------
# PCA on the seeded 70% split, everyone projected, PC scores weighted by variance explained, k-means.
fit_subtype_pipeline <- function(pgs_matrix, k, seed, train_frac = TRAIN_FRAC, nstart = KMEANS_NSTART) {
    set.seed(seed)
    train_idx <- sample(nrow(pgs_matrix), floor(train_frac * nrow(pgs_matrix)))
    pca <- prcomp(pgs_matrix[train_idx, , drop = FALSE], scale. = TRUE)
    variance_explained <- summary(pca)$importance[2, ]
    pc_scores <- predict(pca, pgs_matrix)
    weighted_pcs <- sweep(pc_scores, 2, variance_explained, FUN = '*')
    set.seed(seed)
    kmeans_fit <- kmeans(weighted_pcs, centers = k, nstart = nstart, iter.max = 50)
    list(pca = pca, variance_explained = variance_explained, pc_scores = pc_scores,
         weighted_pcs = weighted_pcs, kmeans = kmeans_fit, train_idx = train_idx, pgs_cols = colnames(pgs_matrix))
}

assign_to_nearest_centroid <- function(fit, pgs_matrix) {
    projected <- sweep(predict(fit$pca, newdata = pgs_matrix), 2, fit$variance_explained, FUN = '*')
    centroids <- fit$kmeans$centers
    squared_distance <- sapply(seq_len(nrow(centroids)), function(i) rowSums(sweep(projected, 2, centroids[i, ])^2))
    apply(squared_distance, 1, which.min)
}

jaccard_similarity <- function(set_a, set_b) {
    if (length(set_a) == 0 && length(set_b) == 0) return(1)
    length(intersect(set_a, set_b)) / length(union(set_a, set_b))
}

# Recode `clusters` to the labels of `reference`, one-to-one, greedily by largest overlap.
match_clusters_one_to_one <- function(clusters, reference) {
    contingency <- table(clusters, reference)
    mapping <- setNames(rep(NA_integer_, nrow(contingency)), rownames(contingency))
    for (step in seq_len(min(dim(contingency)))) {
        largest <- which(contingency == max(contingency), arr.ind = TRUE)[1, ]
        mapping[rownames(contingency)[largest[1]]] <- as.integer(colnames(contingency)[largest[2]])
        contingency[largest[1], ] <- -1L
        contingency[, largest[2]] <- -1L
    }
    unname(mapping[as.character(clusters)])
}

# ---- The published subtypes ---------------------------------------------------------------
subtype_pgs_matrix <- function() memo('subtype_pgs_matrix', {
    cohort <- load_cohort()
    m <- as.matrix(cohort[, pgs_variable_ids()])
    rownames(m) <- cohort[[ID_COL]]
    m
})

fit_published_subtypes <- function() memo('published_subtypes', {
    cohort <- load_cohort()
    fit <- fit_subtype_pipeline(subtype_pgs_matrix(), k = N_SUBTYPES, seed = SPLIT_SEED)
    if (!identical(as.integer(fit$train_idx), as.integer(fit_pgs_pca()$train_idx)))
        stop('fit_published_subtypes: training split differs from fit_pgs_pca()')
    raw <- fit$kmeans$cluster
    gpath_means <- tapply(cohort$gpath, raw, mean, na.rm = TRUE)
    renumber <- setNames(as.integer(rank(-gpath_means, ties.method = 'first')), names(gpath_means))
    subtype <- unname(renumber[as.character(raw)])
    say('Subtypes 1..6 (mean gpath {paste(sprintf("%.3f", sort(gpath_means, decreasing = TRUE)), collapse = "/")}): ',
        'sizes {paste(table(factor(subtype, levels = 1:6)), collapse = "/")}')
    list(subtype = subtype, raw_cluster = unname(raw), fit = fit, colours = subtype_colours())
})

standardized_phenotypes <- function(ids = SUBTYPE_PHENOTYPE_IDS) {
    cohort <- load_cohort()
    out <- as.data.frame(lapply(ids, function(v) as.numeric(scale(cohort[[v]]))))
    names(out) <- ids
    out
}

format_p <- function(p) ifelse(p < 0.001, formatC(p, digits = 1, format = 'e'),
                               formatC(signif(p, 2), digits = 2, format = 'fg', flag = '#'))
p_stars <- function(p) ifelse(p < 0.001, '***', ifelse(p < 0.01, '**', ifelse(p < 0.05, '*', '')))
black_text_theme <- function(size) {
    theme(text = element_text(colour = 'black', size = size), axis.text = element_text(colour = 'black', size = size),
          axis.title = element_text(colour = 'black', size = size), legend.text = element_text(colour = 'black', size = size),
          legend.title = element_text(colour = 'black', size = size), plot.title = element_text(colour = 'black', size = size))
}

# ---- Subtype biplot, shared by Fig 5f and Fig S13 --------------------------------------------
loading_arrows <- function(fit, scale, wrap_width = 26, nudge = 0.4) {
    rotation <- fit$pca$rotation[, 1:2]
    names_ <- mapper_pgs(rownames(rotation))
    wrap <- function(text) paste(strwrap(text, width = wrap_width), collapse = '\n')
    loadings <- data.frame(PC1 = rotation[, 1] * scale, PC2 = rotation[, 2] * scale,
                           label = ifelse(grepl(UNLABELLED_PGS_PATTERN, names_, ignore.case = TRUE), NA_character_, vapply(names_, wrap, character(1))),
                           row.names = rownames(rotation), stringsAsFactors = FALSE)
    labelled <- loadings[!is.na(loadings$label), ]
    len <- sqrt(labelled$PC1^2 + labelled$PC2^2)
    labelled$nudge_x <- labelled$PC1 / len * nudge
    labelled$nudge_y <- labelled$PC2 / len * nudge
    list(all = loadings, labelled = labelled)
}

subtype_biplot <- function(fit, subtypes, colours, title = NULL, legend_title = 'Subtype',
                           legend_labels = function(counts) paste0(1:6, ' (n=', counts, ')'),
                           axis_prefix = 'PC', font_size = HOUSE_FONT_PT, font_family = HOUSE_FONT_FAMILY,
                           point_size = 0.7, arrow_scale = NULL, legend_ncol = 3, limit_factor = 1.32,
                           base_theme = NULL, arrow_linewidth = 0.3, arrow_length = unit(0.13, 'cm'),
                           legend_byrow = TRUE, wrap_width = 26, nudge = 0.4, repel_force = 3, repel_seed = 7) {
    if (is.null(base_theme)) base_theme <- theme_bw(base_size = font_size, base_family = font_family) + black_text_theme(font_size)
    coordinates <- data.frame(PC1 = fit$pc_scores[, 1], PC2 = fit$pc_scores[, 2], subtype = factor(subtypes, levels = 1:6))
    if (is.null(arrow_scale)) arrow_scale <- max(abs(as.matrix(coordinates[, 1:2]))) / max(abs(fit$pca$rotation[, 1:2]))
    arrows <- loading_arrows(fit, arrow_scale, wrap_width = wrap_width, nudge = nudge)
    labelled <- arrows$labelled
    pct <- round(100 * fit$variance_explained[1:2], 1)
    lim_x <- max(abs(c(coordinates$PC1, arrows$all$PC1))) * limit_factor
    lim_y <- max(abs(c(coordinates$PC2, arrows$all$PC2))) * limit_factor
    counts <- as.integer(table(coordinates$subtype))

    ggplot(coordinates, aes(PC1, PC2, colour = subtype)) +
        geom_point(size = point_size, alpha = 0.7) +
        geom_segment(data = labelled, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                     arrow = arrow(length = arrow_length), linewidth = arrow_linewidth, colour = 'black', inherit.aes = FALSE) +
        geom_text_repel(data = labelled, aes(x = PC1, y = PC2, label = label), inherit.aes = FALSE,
                        colour = 'black', family = font_family, size = font_size / .pt,
                        lineheight = 0.85, nudge_x = labelled$nudge_x, nudge_y = labelled$nudge_y,
                        box.padding = 0.45, point.padding = 0.25, force = repel_force, max.overlaps = Inf,
                        min.segment.length = 0.5, segment.colour = 'grey55', segment.size = 0.2, seed = repel_seed) +
        scale_colour_manual(values = colours, name = legend_title, drop = FALSE, labels = legend_labels(counts)) +
        guides(colour = guide_legend(override.aes = list(size = 1.6, alpha = 1), ncol = legend_ncol, byrow = legend_byrow)) +
        labs(title = title, x = glue('{axis_prefix}1 ({pct[1]}%)'), y = glue('{axis_prefix}2 ({pct[2]}%)')) +
        coord_cartesian(xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y)) +
        base_theme +
        theme(legend.key.size = unit(0.32, 'cm'), legend.margin = margin(1, 1, 1, 1), legend.position = 'bottom')
}

# ---- Figure 5f (house style; the legend is shipped as its own PDF) ---------------------------
build_fig5f <- function(with_legend = FALSE, sub = fit_published_subtypes()) {
    plot <- subtype_biplot(sub$fit, sub$subtype, sub$colours,
                           legend_title = 'PGS-based AD polygenic subtypes',
                           legend_labels = function(counts) paste0(1:6, ' (n=', counts, ' individuals)'),
                           axis_prefix = 'PGS PC', point_size = 0.6, arrow_scale = LOADING_ARROW_SCALE, legend_ncol = 2,
                           legend_byrow = FALSE, limit_factor = 1.05, base_theme = theme_house(),
                           arrow_linewidth = HOUSE_LINE_MM, arrow_length = unit(1.5, 'mm'),
                           wrap_width = 20, nudge = 0.5, repel_force = 6, repel_seed = 7) +
        guides(colour = guide_legend(override.aes = list(size = 1.2, alpha = 1), ncol = 2, byrow = FALSE, title.position = 'top'))
    if (!with_legend) return(plot + theme(legend.position = 'none'))
    plot + theme(legend.position = 'bottom', legend.key.size = unit(HOUSE_FONT_PT, 'pt'),
                 legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(2, 'pt'), legend.key.spacing.y = unit(1, 'pt'))
}

# ---- Figure 5g / S12 / S14: plaque violins, VIOLIN_SUBTYPES[1] vs [2] ------------------------
build_plaque_violin <- function(keep = NULL, title = NULL, font_size = 11, base_theme = NULL,
                                linewidth = 0.5, outlier_size = 0.8, bracket_linewidth = 0.6, sub = fit_published_subtypes(),
                                ylim = c(-2.1, 4.5), bracket_y = 4.25) {
    if (is.null(base_theme)) base_theme <- theme_bw(base_size = font_size) + black_text_theme(font_size)
    z <- standardized_phenotypes(names(PLAQUE_IDS))
    if (is.null(keep)) keep <- rep(TRUE, length(sub$subtype))
    in_group <- lapply(VIOLIN_SUBTYPES, function(s) keep & sub$subtype == s)
    names(in_group) <- as.character(VIOLIN_SUBTYPES)
    n <- vapply(in_group, sum, integer(1))
    p <- vapply(names(PLAQUE_IDS), function(v) wilcox.test(z[[v]][in_group[[1]]], z[[v]][in_group[[2]]])$p.value, numeric(1))

    long <- do.call(rbind, lapply(names(PLAQUE_IDS), function(v) {
        do.call(rbind, lapply(as.character(VIOLIN_SUBTYPES), function(s) {
            data.frame(type = PLAQUE_IDS[[v]], subtype = s, value = z[[v]][in_group[[s]]])
        }))
    }))
    long <- long[!is.na(long$value), ]
    long <- long[long$value >= -2 & long$value <= 4, ]
    long$type <- factor(long$type, levels = unname(PLAQUE_IDS))
    long$x <- factor(paste0('Subtype ', long$subtype, '\n(n=', n[long$subtype], ')'),
                     levels = paste0('Subtype ', VIOLIN_SUBTYPES, '\n(n=', n, ')'))
    annotations <- data.frame(type = factor(unname(PLAQUE_IDS), levels = unname(PLAQUE_IDS)),
                              x = 1, xend = 2, x_mid = 1.5, y_line = bracket_y, label = paste0('p=', format_p(p)), stars = p_stars(p))

    plot <- ggplot(long, aes(x = x, y = value, fill = subtype)) +
        geom_violin(trim = FALSE, alpha = 0.7, width = 0.9, linewidth = linewidth) +
        geom_boxplot(width = 0.2, alpha = 0.8, outlier.size = outlier_size, linewidth = linewidth, outlier.stroke = linewidth) +
        geom_segment(data = annotations, aes(x = x, xend = xend, y = y_line, yend = y_line),
                     inherit.aes = FALSE, linewidth = bracket_linewidth, colour = 'black') +
        geom_text(data = annotations, aes(x = x_mid, y = y_line - 0.2, label = label), inherit.aes = FALSE,
                  size = font_size / .pt, vjust = 1, colour = 'black', family = base_theme$text$family) +
        geom_text(data = annotations, aes(x = x_mid, y = y_line + 0.05, label = stars), inherit.aes = FALSE,
                  size = font_size / .pt, vjust = 0, colour = 'black', family = base_theme$text$family) +
        facet_wrap(~type, nrow = 1, strip.position = 'bottom') +
        scale_fill_manual(values = sub$colours[as.character(VIOLIN_SUBTYPES)], guide = 'none') +
        scale_y_continuous(breaks = seq(-2, 4, 2)) +
        coord_cartesian(ylim = ylim, clip = 'off') +
        labs(x = NULL, y = 'Standardized multiregion plaque burden (z-score)', title = title) +
        base_theme +
        theme(strip.placement = 'outside', strip.background = element_blank(),
              strip.text = element_text(colour = 'black', size = font_size),
              panel.spacing = unit(0, 'pt'), plot.title = element_text(hjust = 0.5))
    list(plot = plot, p = p, n = n)
}

# ---- Figure 5h: APOE allelotype per subtype --------------------------------------------------
build_fig5h <- function(sub = fit_published_subtypes()) {
    cohort <- load_cohort()
    apoe <- factor(as.character(cohort$apoe_genotype), levels = names(APOE_LEVELS))
    if (anyNA(apoe)) stop('build_fig5h: apoe_genotype values outside APOE_LEVELS')
    counts <- as.data.frame(table(subtype = factor(sub$subtype, levels = 1:6), apoe = apoe))
    palette <- setNames(rev(RColorBrewer::brewer.pal(6, 'RdBu')), names(APOE_LEVELS))
    ggplot(counts, aes(x = subtype, y = Freq, fill = apoe)) +
        geom_col(position = position_stack(reverse = TRUE), width = 0.85) +
        scale_fill_manual(values = palette, labels = APOE_LEVELS, name = expression(italic(APOE)~'allelotype')) +
        guides(fill = guide_legend(nrow = 1, title.position = 'top', title.hjust = 0.5)) +
        labs(x = 'AD polygenic subtype', y = 'Number of ROSMAP individuals') +
        theme_house() +
        theme(legend.position = 'bottom', legend.key.size = unit(HOUSE_FONT_PT, 'pt'),
              legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, 'pt'),
              legend.spacing.x = unit(1, 'pt'), legend.key.spacing.x = unit(1, 'pt'),
              legend.text = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 1, unit = 'pt')))
}

make_figure5_subtypes <- function() {
    say('Figure 5f-h: AD polygenic subtypes')
    sub <- fit_published_subtypes()

    f_path <- fig_path('Figure5', 'Fig5f_subtype_biplot_PC1_PC2.pdf')
    open_panel_pdf(f_path, 'f'); print(build_fig5f()); dev.off()

    legend <- cowplot::get_plot_component(build_fig5f(with_legend = TRUE), 'guide-box-bottom')
    legend_w <- convertWidth(sum(legend$widths), 'in', valueOnly = TRUE) + 0.04
    legend_h <- convertHeight(sum(legend$heights), 'in', valueOnly = TRUE) + 0.04
    legend_path <- fig_path('Figure5', 'Fig5f_legend_subtypes.pdf')
    open_pdf(legend_path, legend_w, legend_h); grid.newpage(); grid.draw(legend); dev.off()

    violin <- build_plaque_violin(font_size = HOUSE_FONT_PT, base_theme = theme_house(), linewidth = HOUSE_LINE_MM,
                                  outlier_size = 0.3, bracket_linewidth = HOUSE_LINE_MM)
    violin$plot <- violin$plot + labs(y = 'Standardized multiregion\nplaque burden (z-score)')
    g_path <- fig_path('Figure5', glue('Fig5g_plaque_violin_subtype{VIOLIN_SUBTYPES[1]}_vs_{VIOLIN_SUBTYPES[2]}.pdf'))
    open_panel_pdf(g_path, 'g'); print(violin$plot); dev.off()
    say('Fig 5g Wilcoxon p: diffuse {format_p(violin$p["plaq_d"])}, neuritic {format_p(violin$p["plaq_n"])}')

    h_path <- fig_path('Figure5', 'Fig5h_apoe_allelotype_by_subtype.pdf')
    open_panel_pdf(h_path, 'h'); print(build_fig5h()); dev.off()

    rbind(
        manifest_row('Figure5', 'f', f_path, note = glue('k-means k={N_SUBTYPES} on variance-weighted PGS PCs, seed {SPLIT_SEED}')),
        manifest_row('Figure5', 'f (legend)', legend_path, note = glue('standalone legend, {round(legend_w, 2)} x {round(legend_h, 2)} in')),
        manifest_row('Figure5', 'g', g_path, note = glue('Wilcoxon p diffuse={format_p(violin$p["plaq_d"])}, neuritic={format_p(violin$p["plaq_n"])}; n={paste(violin$n, collapse = "/")}')),
        manifest_row('Figure5', 'h', h_path))
}

# ---- Figure S11: subtype means of standardized phenotypes ------------------------------------
make_figureS11 <- function() {
    say('Figure S11: subtype phenotype-mean heatmap')
    sub <- fit_published_subtypes()
    z <- standardized_phenotypes(SUBTYPE_PHENOTYPE_IDS)
    means <- as.matrix(sapply(SUBTYPE_PHENOTYPE_IDS, function(v) tapply(z[[v]], factor(sub$subtype, levels = 1:6), mean, na.rm = TRUE)))
    rownames(means) <- as.character(1:6)
    colnames(means) <- mapper_pheno(SUBTYPE_PHENOTYPE_IDS)
    boxed_cols <- mapper_pheno(names(PLAQUE_IDS))
    boxed_rows <- as.character(VIOLIN_SUBTYPES)

    ht <- Heatmap(means, name = 'subtype_mean', col = colorRamp2(c(-0.4, 0, 0.4), c('violet', 'white', 'red')),
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  row_names_side = 'left', column_names_side = 'bottom', column_names_rot = 90,
                  row_title = 'AD polygenic subtype', row_title_side = 'left',
                  column_title = 'Selected observed ROSMAP AD phenotype', column_title_side = 'bottom',
                  row_names_gp = gpar(fontsize = 11), column_names_gp = gpar(fontsize = 11), column_names_max_height = unit(9, 'cm'),
                  row_title_gp = gpar(fontsize = 12), column_title_gp = gpar(fontsize = 12),
                  heatmap_legend_param = list(title = 'Subtype average of standardized\nobserved AD phenotype (z-score)',
                                              at = c(-0.4, 0, 0.4), labels = c('-0.4', '0', '0.4'),
                                              title_gp = gpar(fontsize = 11), labels_gp = gpar(fontsize = 10)),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                      if (rownames(means)[i] %in% boxed_rows && colnames(means)[j] %in% boxed_cols)
                          grid.rect(x, y, width * 0.72, height * 0.72, gp = gpar(col = '#2C7BB6', fill = NA, lwd = 2.5))
                  })
    path <- fig_path('FigureS11', 'FigS11_subtype_phenotype_means_heatmap.pdf')
    open_pdf(path, width = 6.5, height = 8)
    draw(ht, heatmap_legend_side = 'right')
    dev.off()
    rbind(manifest_row('FigureS11', 'a', path, note = 'means of whole-cohort z-scores per subtype'),
          make_figureS11_pairwise())
}

# ---- Figure S11b-d: all-by-all subtype Wilcoxon tests, FDR over all phenotypes x pairs -------
# Wilcoxon over every subtype pair for each phenotype; FDR across all phenotypes x pairs.
pairwise_subtype_tests <- function(ids = PAIRWISE_PHENOTYPE_IDS, sub = fit_published_subtypes()) {
    z <- standardized_phenotypes(ids)
    pairs <- t(combn(1:6, 2))
    tests <- do.call(rbind, lapply(ids, function(v) {
        stats <- t(apply(pairs, 1, function(ab) {
            x <- z[[v]][sub$subtype == ab[1]]; y <- z[[v]][sub$subtype == ab[2]]
            c(mean_diff = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE), p = wilcox.test(x, y)$p.value)
        }))
        data.frame(phenotype = v, a = pairs[, 1], b = pairs[, 2], stats)
    }))
    tests$fdr <- p.adjust(tests$p, method = 'fdr')
    tests
}

# Half-heatmaps laid on their side: pair (a, b) is a diamond centred at x = (a + b) / 2, height (b - a) / 2,
# coloured by mean z_a - mean z_b, asterisks for FDR. The colour bar is cut at the most negative difference (to 0.1).
pairwise_triangles_plot <- function(tests) {
    ids <- factor(tests$phenotype, levels = PAIRWISE_PHENOTYPE_IDS)
    legend_min <- min(0, floor(min(tests$mean_diff) * 10) / 10)
    ramp <- colorRamp2(c(-PAIRWISE_DIFF_CAP, 0, PAIRWISE_DIFF_CAP), c('violet', 'white', '#DC0000FF'))
    x <- (tests$a + tests$b) / 2; y <- (tests$b - tests$a) / 2
    diamonds <- data.frame(phenotype = rep(ids, each = 4), id = rep(seq_len(nrow(tests)), each = 4),
                           mean_diff = rep(tests$mean_diff, each = 4),
                           px = c(rbind(x - 0.5, x, x + 0.5, x)), py = c(rbind(y, y + 0.5, y, y - 0.5)))
    stars <- data.frame(phenotype = ids, x = x, y = y - 0.1, label = p_stars(tests$fdr))
    zigzag <- rev(seq(1, 6, by = 0.5))[-1]
    outline <- data.frame(px = c(1, 3.5, 6, zigzag), py = c(0.5, 3, 0.5, ifelse(zigzag %% 1 == 0, 0.5, 0)))
    ticks <- expand.grid(x = 1:6, phenotype = factor(PAIRWISE_PHENOTYPE_IDS, levels = PAIRWISE_PHENOTYPE_IDS))

    fill_scale <- function(...) scale_fill_gradientn(colours = ramp(seq(legend_min, PAIRWISE_DIFF_CAP, length.out = 60)),
                                                     limits = c(legend_min, PAIRWISE_DIFF_CAP), guide = 'none', ...)
    triangles <- ggplot() +
        geom_polygon(data = diamonds, aes(px, py, group = id, fill = mean_diff), colour = NA) +
        geom_polygon(data = outline, aes(px, py), fill = NA, colour = 'grey75', linewidth = HOUSE_LINE_MM / 2) +
        geom_text(data = stars, aes(x, y, label = label), size = FIG5_TEXT_SIZE, family = HOUSE_FONT_FAMILY) +
        geom_text(data = ticks, aes(x, -0.12, label = x), size = FIG5_TEXT_SIZE, family = HOUSE_FONT_FAMILY, vjust = 1) +
        fill_scale() +
        facet_wrap(~phenotype, ncol = 1, labeller = as_labeller(setNames(mapper_pheno(PAIRWISE_PHENOTYPE_IDS), PAIRWISE_PHENOTYPE_IDS))) +
        coord_fixed(clip = 'off') + scale_y_continuous(limits = c(-0.45, 3)) +
        labs(x = 'AD polygenic subtype', y = NULL,
             caption = glue('FDR: * < 0.05, ** < 0.01, *** < 0.001 (Wilcoxon, FDR over {nrow(tests)} tests)')) +
        theme_void(base_size = HOUSE_FONT_PT, base_family = HOUSE_FONT_FAMILY) +
        theme(text = element_text(size = HOUSE_FONT_PT, colour = HOUSE_TEXT_COLOUR),
              strip.text = element_text(size = HOUSE_FONT_PT, margin = margin(b = 3)),
              axis.title.x = element_text(size = HOUSE_FONT_PT, margin = margin(t = 3)),
              plot.caption = element_text(size = HOUSE_FONT_PT, hjust = 0), plot.margin = margin(6, 111, 6, 6), panel.spacing = unit(9, 'pt'))

    # Colour bar drawn by hand (S11a-sized): ends labelled on the right, 0 on the left. Bar spans x 0..0.5.
    step <- (PAIRWISE_DIFF_CAP - legend_min) / 100
    bar <- data.frame(y = seq(legend_min + step / 2, PAIRWISE_DIFF_CAP - step / 2, by = step))
    right_ticks <- data.frame(y = c(legend_min, PAIRWISE_DIFF_CAP))
    colour_bar <- ggplot() +
        geom_tile(data = bar, aes(0.25, y, fill = y), width = 0.5, height = step) + fill_scale() +
        geom_segment(data = right_ticks, aes(x = 0.5, xend = 0.62, y = y, yend = y), linewidth = HOUSE_LINE_MM) +
        geom_text(data = right_ticks, aes(0.7, y, label = y), hjust = 0, size = FIG5_TEXT_SIZE, family = HOUSE_FONT_FAMILY) +
        annotate('segment', x = -0.12, xend = 0, y = 0, yend = 0, linewidth = HOUSE_LINE_MM) +
        annotate('text', x = -0.2, y = 0, label = '0', hjust = 1, size = FIG5_TEXT_SIZE, family = HOUSE_FONT_FAMILY) +
        scale_x_continuous(limits = c(-0.6, 5), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
        coord_cartesian(clip = 'off') +
        labs(title = 'Difference in subtype average of\nstandardized phenotype\n(left - right subtype, z-score)') +
        theme_void(base_size = HOUSE_FONT_PT, base_family = HOUSE_FONT_FAMILY) +
        theme(plot.title = element_text(size = HOUSE_FONT_PT, colour = HOUSE_TEXT_COLOUR, hjust = 0, margin = margin(b = 4)),
              plot.title.position = 'plot')
    triangles + inset_element(colour_bar, left = 0.64, right = 1, bottom = 0.42, top = 0.42 + PAIRWISE_BAR_FRACTION,
                              align_to = 'full', clip = FALSE)
}

make_figureS11_pairwise <- function() {
    say('Figure S11b-d: all-by-all subtype comparisons ({paste(PAIRWISE_PHENOTYPE_IDS, collapse = ", ")})')
    tests <- pairwise_subtype_tests()
    path <- fig_path('FigureS11', 'FigS11bcd_pairwise_subtypes.pdf')
    ggsave(path, pairwise_triangles_plot(tests), width = S11_PAIRWISE_IN[1], height = S11_PAIRWISE_IN[2], device = open_pdf)
    n_sig <- vapply(PAIRWISE_PHENOTYPE_IDS, function(v) sum(tests$phenotype == v & tests$fdr < 0.05), integer(1))
    say('Fig S11b-d: pairs with FDR < 0.05 (of 15): {paste(names(n_sig), n_sig, sep = " ", collapse = ", ")}')
    manifest_row('FigureS11', 'b-d', path, note = glue('pairwise Wilcoxon, FDR over {nrow(tests)} tests; pairs FDR < 0.05 of 15: ',
                                                      '{paste(names(n_sig), n_sig, sep = "=", collapse = ", ")}'))
}

# ---- Figure S12: Fig 5g within ε3/ε3 individuals --------------------------------------------
make_figureS12 <- function() {
    say('Figure S12: plaque violins within ε3/ε3 individuals')
    violin <- build_plaque_violin(keep = load_cohort()$apoe_genotype == 0, title = 'ε3/ε3 individuals')
    path <- fig_path('FigureS12', 'FigS12_plaque_violin_e3e3.pdf')
    ggsave(path, violin$plot, width = 5, height = 4.5, device = open_pdf)
    say('Fig S12 Wilcoxon p: diffuse {format_p(violin$p["plaq_d"])}, neuritic {format_p(violin$p["plaq_n"])}; n={paste(violin$n, collapse = "/")}')
    manifest_row('FigureS12', 'violin', path, note = glue('apoe_genotype == 0; Wilcoxon p diffuse={format_p(violin$p["plaq_d"])}, ',
                                                         'neuritic={format_p(violin$p["plaq_n"])}; n={paste(violin$n, collapse = "/")}'))
}

# ---- Bootstrap resamples (Fig S13, Table S11) -----------------------------------------------
# Refit on a bootstrap resample, relabel everyone by nearest centroid, match to the published numbering.
resample_subtypes <- function(iteration) {
    pgs_matrix <- subtype_pgs_matrix()
    reference <- fit_published_subtypes()$subtype
    seed <- SPLIT_SEED + iteration
    set.seed(seed)
    resample_idx <- sample(seq_len(nrow(pgs_matrix)), size = nrow(pgs_matrix), replace = TRUE)
    fit <- fit_subtype_pipeline(pgs_matrix[resample_idx, , drop = FALSE], k = N_SUBTYPES, seed = seed)
    match_clusters_one_to_one(assign_to_nearest_centroid(fit, pgs_matrix), reference)
}

subtype_jaccards <- function(resampled, reference = fit_published_subtypes()$subtype) {
    sapply(1:6, function(s) jaccard_similarity(which(reference == s), which(resampled == s)))
}

build_figS13 <- function(sub = fit_published_subtypes()) {
    panels <- c(
        list(subtype_biplot(sub$fit, sub$subtype, sub$colours, title = 'Primary analysis')),
        lapply(seq_along(S13_EXAMPLE_ITERATIONS), function(i) {
            resampled <- resample_subtypes(S13_EXAMPLE_ITERATIONS[i])
            jaccard <- mean(subtype_jaccards(resampled, sub$subtype))
            say('S13 resample {i} (iteration {S13_EXAMPLE_ITERATIONS[i]}): mean Jaccard {round(jaccard, 4)}')
            subtype_biplot(sub$fit, resampled, sub$colours, title = glue('Resample {i} (mean Jaccard = {sprintf("%.2f", jaccard)})'))
        }))
    wrap_plots(panels, ncol = 2)
}

# Table S11: mean and SD over the resamples of each subtype's Jaccard stability.
bootstrap_stability_table <- function(n_bootstrap = N_BOOTSTRAP, reference = fit_published_subtypes()$subtype) {
    jaccard_matrix <- t(sapply(seq_len(n_bootstrap), function(b) subtype_jaccards(resample_subtypes(b), reference)))
    data.frame(`AD Polygenic Subtype` = 1:6,
               `Mean Bootstrap Jaccard Stability` = colMeans(jaccard_matrix),
               `Jaccard Stability Bootstrap Standard Deviation` = apply(jaccard_matrix, 2, sd), check.names = FALSE)
}

make_figureS13_and_tables <- function(n_bootstrap = N_BOOTSTRAP) {
    say('Figure S13: bootstrap resample examples')
    s13_path <- fig_path('FigureS13', 'FigS13_bootstrap_resample_examples.pdf')
    ggsave(s13_path, build_figS13(), width = 7.2, height = 7.6, device = open_pdf)

    say('Table S11: bootstrap stability, {n_bootstrap} resamples')
    s11_path <- table_path('TableS11_bootstrap_stability_k6.tsv')
    write.table(bootstrap_stability_table(n_bootstrap), s11_path, sep = '\t', quote = FALSE, row.names = FALSE)

    rbind(
        manifest_row('FigureS13', 'a-d', s13_path, note = glue('primary fit + bootstrap iterations {paste(S13_EXAMPLE_ITERATIONS, collapse = ", ")}')),
        manifest_row('TableS11', 'bootstrap stability', s11_path,
                     note = if (n_bootstrap == N_BOOTSTRAP) glue('{n_bootstrap} resamples') else glue('QUICK run: {n_bootstrap} resamples (the paper used {N_BOOTSTRAP})')))
}

# ============================ PART 3: Figure S14 ======================================
# k-means on the training rows only; test rows by nearest centroid; labels matched on the training rows.
fit_subtypes_holdout <- function(seed = SPLIT_SEED) memo(glue('holdout_subtypes_{seed}'), {
    published <- fit_published_subtypes()
    pca_fit <- fit_pgs_pca()
    train <- pca_fit$train_idx; test <- pca_fit$test_idx
    if (!identical(as.integer(train), as.integer(published$fit$train_idx)))
        stop('fit_subtypes_holdout: fit_pgs_pca() training split differs from the published subtype split')

    weighted <- sweep(pca_fit$scores_all, 2, pca_fit$variance_explained, FUN = '*')
    set.seed(seed)
    km <- kmeans(weighted[train, , drop = FALSE], centers = N_SUBTYPES, nstart = KMEANS_NSTART, iter.max = 50)
    squared_distance <- sapply(seq_len(N_SUBTYPES), function(i) rowSums(sweep(weighted, 2, km$centers[i, ])^2))
    raw <- apply(squared_distance, 1, which.min)
    if (!all(raw[train] == km$cluster)) stop('fit_subtypes_holdout: nearest-centroid rule does not reproduce the training clusters')

    mapping <- match_clusters_one_to_one(km$cluster, published$subtype[train])
    lookup <- setNames(mapping, km$cluster)[!duplicated(km$cluster)]
    subtype <- unname(as.integer(lookup[as.character(raw)]))
    if (!identical(sort(unique(subtype)), 1:6)) stop('fit_subtypes_holdout: matching did not yield subtypes 1..6')

    agreement <- c(train = mean(subtype[train] == published$subtype[train]), test = mean(subtype[test] == published$subtype[test]))
    say('Hold-out subtypes: train sizes {paste(table(factor(subtype[train], levels = 1:6)), collapse = "/")}, ',
        'test sizes {paste(table(factor(subtype[test], levels = 1:6)), collapse = "/")}; agreement with published labels ',
        '{sprintf("%.1f", 100 * agreement["train"])}% (train) / {sprintf("%.1f", 100 * agreement["test"])}% (test)')
    list(subtype = subtype, raw_cluster = unname(raw), kmeans = km, train_idx = train, test_idx = test,
         agreement = agreement, colours = subtype_colours(), published = published$subtype)
})

make_figureS14 <- function() {
    say('Figure S14: Fig 5g on the held-out test split, subtypes fit on the training split')
    sub <- fit_subtypes_holdout()
    label <- glue('Held-out test split (30%, n={length(sub$test_idx)})')
    violin <- build_plaque_violin(keep = seq_along(sub$subtype) %in% sub$test_idx, title = label, sub = sub,
                                  ylim = c(-2.8, 6), bracket_y = 5.6)
    path <- fig_path('FigureS14', 'FigS14_plaque_violin_test_split.pdf')
    ggsave(path, violin$plot, width = 6, height = 5.5, device = open_pdf)
    say('Fig S14 Wilcoxon p: diffuse {format_p(violin$p["plaq_d"])}, neuritic {format_p(violin$p["plaq_n"])}; n={paste(violin$n, collapse = "/")}')
    manifest_row('FigureS14', 'violin', path, note = glue('{label}; k-means on training rows only; Wilcoxon p diffuse={format_p(violin$p["plaq_d"])}, ',
                                                         'neuritic={format_p(violin$p["plaq_n"])}; n={paste(violin$n, collapse = "/")}'))
}
