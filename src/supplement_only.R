# src/supplement_only.R -- supplementary figures not tied to a main-text figure.
#   make_figureS1()      Fig S1  phenotype heatmap, 36 phenotypes x 1678 individuals
#   make_figureS2()      Fig S2  PCA of the 36 mean-imputed phenotypes
#   make_figureS4()      Fig S4  ancestry PCA against the HGDP / 1000 Genomes reference
#   make_figureS6()      Fig S6  AD PGS transferability to NIA-Reagan AD status
#   make_figureS9_S10()  Fig S9  GREAT enrichment of the ApoB and prostate-cancer PGS models
#                        Fig S10 treemap of the ApoB GO terms by semantic (GloVe) cluster
suppressPackageStartupMessages({
    library(ggplot2)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
    library(ggrepel)
    library(pROC)
    library(scales)
})

AD_PGS_COL <- 'FH1263_SCORE1_AVG'

zscore_na <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
# Population (ddof = 0) z-score, as scipy.stats.zscore in the original notebook.
zscore_pop <- function(x) { mu <- mean(x); (x - mu) / sqrt(mean((x - mu)^2)) }
wrap_label <- function(x, width = 24) vapply(x, function(s) paste(strwrap(s, width), collapse = '\n'), character(1))

# ---- Figure S1: phenotype heatmap ------------------------------------------------------------
make_figureS1 <- function() {
    say('Figure S1: phenotype heatmap')
    phenotypes <- load_phenotypes()
    ids <- stage2_phenotype_ids()
    rownames(phenotypes) <- phenotypes[[ID_COL]]
    m <- as.matrix(phenotypes[, ids])
    colnames(m) <- mapper_pheno(colnames(m))

    z <- apply(m, 2, zscore_na)
    imputed <- apply(z, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x })
    z <- t(z); imputed <- t(imputed)
    n_pheno <- nrow(z); n_ind <- ncol(z)

    ht <- Heatmap(z, width = n_ind * unit(0.24, 'mm'), height = n_pheno * unit(12, 'mm'),
                  column_dend_height = unit(2, 'cm'), row_dend_width = unit(2, 'cm'),
                  col = colorRamp2(c(-2, 0, 2), c('violet', 'white', 'red')), na_col = 'grey',
                  column_title = glue('{n_ind} ROSMAP Individuals'), column_title_gp = gpar(fontsize = 22),
                  row_title = glue('{n_pheno} observed ROSMAP AD phenotypes'), row_title_gp = gpar(fontsize = 22),
                  cluster_rows = TRUE, cluster_columns = TRUE,
                  clustering_distance_rows = dist(imputed), clustering_distance_columns = dist(t(imputed)),
                  show_row_names = TRUE, show_column_names = FALSE, row_names_gp = gpar(fontsize = 22),
                  row_names_max_width = max_text_width(colnames(m), gp = gpar(fontsize = 22)),
                  heatmap_legend_param = list(title = 'Standardized ROSMAP AD phenotype (z-score)', legend_direction = 'horizontal',
                                              title_gp = gpar(fontsize = 22), labels_gp = gpar(fontsize = 22)))
    path <- fig_path('FigureS1', glue('FigS1_phenotype_heatmap_{n_pheno}x{n_ind}.pdf'))
    pdf(path, width = 27, height = 21.5)
    draw(ht, heatmap_legend_side = 'top', padding = unit(c(2, 2, 2, 2), 'cm'))
    dev.off()
    manifest_row('FigureS1', 'S1', path, note = glue('{n_pheno} Stage-2 phenotypes x {n_ind} individuals, z-scored per phenotype; grey = missing'))
}

# ---- Figure S2: PCA of the 36 standardized phenotypes -----------------------------------------
make_figureS2 <- function() {
    say('Figure S2: phenotype PCA')
    full_input <- mean_impute(load_cohort())
    pheno_ids <- stage2_phenotype_ids()
    pc <- prcomp(full_input[, pheno_ids], scale. = TRUE)
    var_prop <- summary(pc)$importance[2, ]

    # (a) variance explained, top 12 PCs
    n_show <- 12
    df_var <- data.frame(PC = factor(paste0('PC', 1:n_show), levels = paste0('PC', 1:n_show)), Variance_Explained = var_prop[1:n_show])
    p_var <- ggplot(df_var, aes(x = PC, y = Variance_Explained)) +
        geom_col(fill = '#2B6CB0', colour = 'blue') +
        geom_text(aes(label = percent(Variance_Explained, accuracy = 0.1)), vjust = -0.5, size = 5) +
        labs(x = glue('ROSMAP phenotype principal components (PCs)\n(top {n_show} of {length(pheno_ids)})'),
             y = 'Proportion of variance explained') +
        theme_bw(base_size = 18)

    # (b) PC1 vs PC2, individuals coloured by AD PGS z-score, arrows for seven highlighted phenotypes
    scores <- as.data.frame(pc$x[, 1:2]); colnames(scores) <- c('PC1', 'PC2')
    scores$AD_PGS <- as.numeric(scale(full_input[[AD_PGS_COL]]))
    loadings <- as.data.frame(pc$rotation[, 1:2]); colnames(loadings) <- c('PC1', 'PC2')
    highlight <- c('amyloid', 'nft', 'plaq_n', 'plaq_d', 'tangles', 'gpath', 'cogn_global_lv')
    loadings[!rownames(loadings) %in% highlight, c('PC1', 'PC2')] <- 0
    loadings$var <- ifelse(rownames(loadings) %in% highlight, wrap_label(mapper_pheno(rownames(loadings))), '')
    scale_factor <- 35
    x_lim <- max(abs(c(scores$PC1, loadings$PC1 * scale_factor)))
    y_lim <- max(abs(c(scores$PC2, loadings$PC2 * scale_factor)))
    pct <- function(i) percent(var_prop[i], accuracy = 0.1)
    p_biplot <- ggplot(scores, aes(x = PC1, y = PC2, colour = AD_PGS)) +
        geom_point(size = 1) +
        geom_segment(data = loadings, aes(x = 0, y = 0, xend = PC1 * scale_factor, yend = PC2 * scale_factor),
                     arrow = arrow(length = unit(0.2, 'cm')), colour = 'red') +
        geom_text_repel(data = loadings, aes(x = PC1 * scale_factor, y = PC2 * scale_factor, label = var),
                        size = 5, colour = 'black', max.overlaps = Inf, seed = SPLIT_SEED) +
        scale_colour_gradient(low = 'white', high = 'blue', name = 'AD PGS\n(z-score)') +
        labs(x = glue('ROSMAP phenotype PC1 ({pct(1)})'), y = glue('ROSMAP phenotype PC2 ({pct(2)})')) +
        coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(-y_lim, y_lim), clip = 'off') +
        theme_bw(base_size = 16)

    pa <- fig_path('FigureS2', 'FigS2a_phenotype_pca_variance.pdf')
    pb <- fig_path('FigureS2', 'FigS2b_phenotype_pca_biplot.pdf')
    ggsave(pa, p_var, width = 8, height = 6)
    ggsave(pb, p_biplot, width = 8, height = 6)
    rbind(manifest_row('FigureS2', 'S2a', pa, note = glue('PCA of {length(pheno_ids)} mean-imputed, standardized phenotypes; PC1-3 = {paste(percent(var_prop[1:3], accuracy = 0.1), collapse = "/")}')),
          manifest_row('FigureS2', 'S2b', pb, note = 'PC1 vs PC2 coloured by AD PGS z-score; loadings x35 for 7 highlighted phenotypes'))
}

# ---- Figure S4: ancestry PCA against the HGDP/1000 Genomes reference --------------------------
# plink2 --score projections (scripts/paths.R entries below); skipped when they are not provided.

ANCESTRY_PCA_FILES <- c(reference = 'ancestry_reference', tgen = 'ancestry_tgen', broad = 'ancestry_broad')
SEABORN_DEEP <- c('#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3', '#937860')
SUPERPOP_LEVELS <- c('EUR', 'EAS', 'AMR', 'SAS', 'AFR')
SUPERPOP_LABELS <- c(EUR = 'European', EAS = 'East Asian', AMR = 'Admixed American', SAS = 'South Asian', AFR = 'African')

read_sscore <- function(key) read.delim(input_path(key), header = TRUE, sep = '\t', check.names = TRUE, stringsAsFactors = FALSE)

build_ancestry_panel <- function(batch_key, batch_label, base_size = 21) {
    ref <- read_sscore(ANCESTRY_PCA_FILES[['reference']])
    batch <- read_sscore(batch_key)
    batch$SuperPop <- batch_label
    batch$Population <- NA_character_
    common <- intersect(colnames(ref), colnames(batch))
    combined <- rbind(ref[, common], batch[, common])
    pc_cols <- grep('^PC[0-9]+_(AVG|SUM)$', colnames(combined), value = TRUE)
    combined[pc_cols] <- lapply(combined[pc_cols], zscore_pop)

    levels_all <- c(SUPERPOP_LEVELS, batch_label)
    labels_all <- c(SUPERPOP_LABELS[SUPERPOP_LEVELS], batch_label)
    combined$SuperPop <- factor(combined$SuperPop, levels = levels_all)
    axis_lab <- function(i) glue('Human Genome Diversity Project/1000 Genomes Project\nprincipal component {i} (z-score)')
    ggplot(combined, aes(x = PC1_AVG, y = PC2_AVG, colour = SuperPop)) +
        geom_point(size = 2, shape = 16) +
        scale_colour_manual(values = setNames(SEABORN_DEEP[seq_along(levels_all)], levels_all), labels = setNames(labels_all, levels_all), name = NULL) +
        labs(x = axis_lab(1), y = axis_lab(2)) +
        guides(colour = guide_legend(override.aes = list(size = 5))) +
        theme_bw(base_size = base_size) +
        theme(panel.grid = element_blank(), panel.border = element_rect(colour = 'black', linewidth = 1),
              legend.position = 'inside', legend.position.inside = c(0.98, 0.02), legend.justification.inside = c(1, 0),
              legend.background = element_rect(colour = 'grey80'), legend.text = element_text(size = base_size))
}

make_figureS4 <- function() {
    say('Figure S4: ancestry PCA')
    pa <- fig_path('FigureS4', 'FigS4a_ancestry_pca_broad.pdf')
    pb <- fig_path('FigureS4', 'FigS4b_ancestry_pca_tgen.pdf')
    missing <- ANCESTRY_PCA_FILES[vapply(ANCESTRY_PCA_FILES, function(k) is.null(optional_path(k)), logical(1))]
    if (length(missing) > 0) {
        note <- glue('ancestry inputs not provided in scripts/paths.R: {paste(missing, collapse = ", ")}')
        return(rbind(manifest_row('FigureS4', 'S4a', pa, 'skipped', note), manifest_row('FigureS4', 'S4b', pb, 'skipped', note)))
    }
    ggsave(pa, build_ancestry_panel(ANCESTRY_PCA_FILES[['broad']], 'Broad Genotyped ROSMAP Individuals'), width = 10, height = 10)
    ggsave(pb, build_ancestry_panel(ANCESTRY_PCA_FILES[['tgen']],  'TGen Genotyped ROSMAP Individuals'),  width = 10, height = 10)
    rbind(manifest_row('FigureS4', 'S4a', pa, note = 'Broad-genotyped ROSMAP + HGDP/1KG reference; PC1/PC2 z-scored over combined rows'),
          manifest_row('FigureS4', 'S4b', pb, note = 'TGen-genotyped ROSMAP + HGDP/1KG reference; same processing'))
}

# ---- Figure S6: AD PGS transferability to NIA-Reagan AD status ---------------------------------
make_figureS6 <- function() {
    say('Figure S6: AD PGS transferability (NIA-Reagan)')
    df <- load_cohort()
    d <- df[!is.na(df$niareagansc) & !is.na(df[[AD_PGS_COL]]), ]
    case_lab <- 'Case\n(NIA-Reagan score <= 2)'
    ctrl_lab <- 'Control\n(NIA-Reagan score > 2)'
    d$status <- factor(ifelse(d$niareagansc <= 2, case_lab, ctrl_lab), levels = c(case_lab, ctrl_lab))
    d$AD_PGS <- as.numeric(scale(d[[AD_PGS_COL]]))
    roc_obj <- roc(cases = d$AD_PGS[d$status == case_lab], controls = d$AD_PGS[d$status == ctrl_lab], quiet = TRUE)
    auroc <- as.numeric(auc(roc_obj))
    n_case <- sum(d$status == case_lab); n_control <- sum(d$status == ctrl_lab)
    say('  n = {nrow(d)}/{nrow(df)} with NIA-Reagan + PGS ({n_case} cases, {n_control} controls); AUROC = {round(auroc, 4)}')

    p_violin <- ggplot(d, aes(x = status, y = AD_PGS)) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width = 0.1) +
        theme_bw(base_size = 18) +
        theme(panel.grid.major.x = element_line(colour = 'grey', linewidth = 0.5), panel.grid.minor.x = element_line(colour = 'grey', linewidth = 0.25)) +
        labs(y = 'AD PGS (z-score)', x = "Alzheimer's disease status based on NIA-Reagan score")
    p_roc <- ggroc(roc_obj, legacy.axes = TRUE, linewidth = 1.2, colour = '#1c61b6') +
        theme_bw(base_size = 18) +
        labs(x = '1 - specificity', y = 'Sensitivity') +
        annotate('text', x = 0.65, y = 0.3, size = 5.5, label = glue('AUROC = {round(auroc, 3)} for the binary\nclassification of AD with AD PGS'))

    pa <- fig_path('FigureS6', 'FigS6a_ad_pgs_by_niareagan_violin.pdf')
    pb <- fig_path('FigureS6', 'FigS6b_ad_pgs_roc.pdf')
    ggsave(pa, p_violin, width = 8, height = 6.5)
    ggsave(pb, p_roc,    width = 8, height = 6.5)
    rbind(manifest_row('FigureS6', 'S6a', pa, note = glue('AD PGS z-scored within the {nrow(d)} individuals with NIA-Reagan; case = score <= 2 (n={n_case}), control > 2 (n={n_control})')),
          manifest_row('FigureS6', 'S6b', pb, note = glue('ROC of AD PGS for NIA-Reagan case vs control; AUROC = {round(auroc, 3)}')))
}

# ---- Figures S9 and S10: GREAT enrichment and semantic clustering of the ApoB GO terms ---------
GREAT_TABLE_KEYS <- c(INI30640 = 'great_apob', cancer1044 = 'great_prostate')

plot_great_enrichment <- function(pgs_id) {
    tb <- read.table(input_path(GREAT_TABLE_KEYS[[pgs_id]]), header = TRUE, sep = '\t')
    tb <- tb[tb$HyperFdrQ < 0.05, ]
    tb$fold_change <- as.numeric(tb$RegionFoldEnrich)
    tb$p_value <- as.numeric(tb$BinomFdrQ)
    tb$delabel <- ifelse(tb$Desc %in% head(tb$Desc, 10), tb$Desc, NA)
    tb$diffexpressed <- 'NO'
    tb$diffexpressed[tb$fold_change > 1.5 & tb$p_value < 0.05] <- 'UP'
    ggplot(data = tb, aes(x = fold_change, y = -log10(p_value), col = diffexpressed, label = delabel)) +
        geom_vline(xintercept = 1.5, col = 'black', linetype = 'dashed') +
        geom_hline(yintercept = -log10(0.05), col = 'black', linetype = 'dashed') +
        theme_bw(base_size = 20) +
        theme(text = element_text(size = 20, color = 'black'), axis.text = element_text(size = 20, color = 'black'),
              axis.title = element_text(size = 20, color = 'black'), legend.text = element_text(size = 20, color = 'black'),
              legend.title = element_text(size = 20, color = 'black'), plot.title = element_text(size = 20, color = 'black', hjust = 0),
              plot.margin = margin(5.5, 5.5, 5.5, 5.5)) +
        geom_point(size = 2) +
        scale_color_manual(values = c('grey', 'blue'), labels = c('Not significant', 'Enriched')) +
        labs(color = 'Severe', x = 'Fold-change', y = expression('-log'[10] * 'q')) +
        geom_text_repel(size = 5, color = 'black') +
        annotate('text', x = max(tb$fold_change, na.rm = TRUE), y = -log10(0.05), label = 'FDR = 0.05', color = 'black', size = 6, hjust = 1) +
        annotate('text', x = 1.5, y = max(-log10(tb$p_value), na.rm = TRUE), label = '1.5 fold-change', angle = 90, color = 'black', size = 6, hjust = 0) +
        scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(1, NA), breaks = c(1, pretty(tb$fold_change[tb$fold_change > 1], n = 5)))
}

GO_THEMES <- c(cholesterol = 'reverse cholesterol transport', lipid = 'lipid localization', catabolism = 'lipid catabolic process',
               transport = 'lipid transport', chylomicron = 'chylomicron assembly', `metabolic processes` = 'cholesterol metabolic process',
               homeostasis = 'lipid homeostasis', lipoprotein = 'plasma lipoprotein particle remodeling')

# The 61 leading ApoB GO terms: HyperFdrQ < 0.05 and BinomFdrQ < 1e-15, in table order.
leading_apob_terms <- function() {
    tb <- read.delim(input_path('great_apob'), comment.char = '#', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    tb[tb$HyperFdrQ < 0.05 & tb$BinomFdrQ < 1e-15, ]
}

# spaCy-like tokenizer: whitespace split, punctuation and infix hyphens as their own tokens.
tokenize_term <- function(s) {
    s <- gsub('([,;:()\\[\\]"])', ' \\1 ', s)
    s <- gsub('(?<=[A-Za-z0-9])-(?=[A-Za-z0-9])', ' - ', s, perl = TRUE)
    toks <- strsplit(trimws(gsub('\\s+', ' ', s)), ' ')[[1]]
    toks[nzchar(toks)]
}

# Mean GloVe vector of each term's tokens (minus of/to/-); streams the vectors file once.
embed_terms <- function(terms, glove_file) {
    tokens <- unique(unlist(lapply(terms, tokenize_term)))
    emb <- list(); con <- file(glove_file, 'r'); on.exit(close(con))
    while (length(lines <- readLines(con, n = 20000, warn = FALSE)) > 0) {
        w <- sub(' .*$', '', lines); hit <- which(w %in% tokens & !(w %in% names(emb)))
        for (i in hit) emb[[w[i]]] <- as.numeric(strsplit(lines[i], ' ')[[1]][-1])
        if (length(emb) == length(tokens)) break
    }
    if (length(emb) == 0) stop('embed_terms: no tokens found in the GloVe file')
    E <- do.call(rbind, emb)
    m <- t(sapply(terms, function(p) {
        tk <- unique(tokenize_term(p)); tk <- tk[tk %in% rownames(E) & !(tk %in% c('of', 'to', '-'))]
        colMeans(E[tk, , drop = FALSE])
    }))
    colnames(m) <- paste0('V', seq_len(ncol(m)))
    m
}

cluster_go_terms <- function(n_clusters = 8) {
    phrase_matrix <- embed_terms(leading_apob_terms()$Desc, input_path('glove_vectors'))
    set.seed(123)
    clusters <- kmeans(phrase_matrix, centers = n_clusters, nstart = 25)$cluster
    anchor_cluster <- clusters[GO_THEMES]
    themes_resolved <- !any(is.na(anchor_cluster)) && !any(duplicated(anchor_cluster)) && length(anchor_cluster) == n_clusters
    if (themes_resolved) {
        labels <- setNames(names(GO_THEMES), anchor_cluster)[as.character(clusters)]
        levels <- names(GO_THEMES)
    } else {
        labels <- paste('Cluster', clusters)
        levels <- paste('Cluster', seq_len(n_clusters))
    }
    list(data = data.frame(phrase = names(clusters), cluster = clusters, cluster_label = factor(labels, levels = levels), size = 1, stringsAsFactors = FALSE),
         phrase_matrix = phrase_matrix, themes_resolved = themes_resolved)
}

make_figureS9_S10 <- function() {
    rows <- list()

    # ---- S9a/b: GREAT enrichment scatters ----------------------------------------------
    scatter_files <- c(INI30640 = 'FigS9a_great_apob_enrichment.pdf', cancer1044 = 'FigS9b_great_prostate_enrichment.pdf')
    panels <- c(INI30640 = 'a', cancer1044 = 'b')
    for (pgs_id in names(GREAT_PGS)) {
        path <- fig_path('FigureS9', scatter_files[[pgs_id]])
        ggsave(path, plot_great_enrichment(pgs_id), width = 12, height = 8, dpi = 300)
        rows[[length(rows) + 1]] <- manifest_row('FigureS9', panels[[pgs_id]], path,
                                                 note = glue('{GREAT_PGS[[pgs_id]]} PGS: GREAT fold enrichment vs -log10 binomial FDR (HyperFdrQ < 0.05)'))
    }

    # ---- S9c / S10 need the GloVe vectors ------------------------------------------------
    s9c_path <- fig_path('FigureS9', 'FigS9c_apob_top25_go_terms.pdf')
    s10_path <- fig_path('FigureS10', 'FigS10_apob_go_treemap.pdf')
    if (is.null(optional_path('glove_vectors'))) {
        note <- 'glove_vectors not provided in scripts/paths.R; semantic clustering skipped'
        rows[[length(rows) + 1]] <- manifest_row('FigureS9', 'c', s9c_path, 'skipped', note)
        rows[[length(rows) + 1]] <- manifest_row('FigureS10', 'treemap', s10_path, 'skipped', note)
        say('Figure S9a-b done; S9c and S10 skipped ({note})')
        return(do.call(rbind, rows))
    }

    clustering <- cluster_go_terms(n_clusters = 8)
    treemap_data <- clustering$data
    cluster_note <- if (clustering$themes_resolved) 'clusters named by anchor term, as published' else
        'WARNING: k-means clusters no longer match the published themes; generic Cluster N labels used'
    # ---- S9c: top-25 ApoB GO terms coloured by cluster -------------------------------------
    df <- leading_apob_terms()[, c('Desc', 'BinomFdrQ')]
    df$neg_log_q <- -log10(df$BinomFdrQ)
    top_25 <- merge(head(df[order(df$neg_log_q, decreasing = TRUE), ], 25), treemap_data[, c('phrase', 'cluster_label')],
                    by.x = 'Desc', by.y = 'phrase', all.x = TRUE)
    p_bar <- ggplot(top_25, aes(x = reorder(Desc, neg_log_q), y = neg_log_q, fill = cluster_label)) +
        geom_bar(stat = 'identity') +
        coord_flip() +
        scale_fill_brewer(palette = 'Set3', name = 'Cluster', drop = FALSE) +
        labs(title = 'Top 25 Gene Ontology Terms by -log10(BinomFdrQ)', x = 'GO Process', y = '-log10(BinomFdrQ)') +
        theme_bw() +
        theme(axis.text.y = element_text(size = 13), axis.text.x = element_text(size = 13))
    ggsave(s9c_path, p_bar, width = 10, height = 8, dpi = 300)
    rows[[length(rows) + 1]] <- manifest_row('FigureS9', 'c', s9c_path, note = glue('top-25 ApoB GO terms by semantic cluster; {cluster_note}'))

    # ---- S10: treemap of all 61 terms by cluster (treemap package, else treemapify) ----------
    palette <- RColorBrewer::brewer.pal(min(nlevels(treemap_data$cluster_label), 11), 'Set3')
    if (requireNamespace('treemap', quietly = TRUE)) {
        pdf(s10_path, width = 12, height = 8)
        treemap::treemap(treemap_data, index = c('cluster_label', 'phrase'), vSize = 'size', type = 'categorical',
                         vColor = 'cluster_label', palette = palette,
                         title = 'Clustered Gene Ontology Terms', title.legend = 'cluster',
                         fontsize.labels = c(0, 12), fontcolor.labels = c('transparent', 'black'), fontface.labels = c(2, 1),
                         align.labels = list(c('center', 'center'), c('left', 'top')), overlap.labels = 0.5, inflate.labels = FALSE)
        invisible(dev.off())
        engine <- 'treemap'
    } else {
        p_tree <- ggplot(treemap_data, aes(area = size, fill = cluster_label, label = phrase, subgroup = cluster_label)) +
            treemapify::geom_treemap(colour = 'black') +
            treemapify::geom_treemap_subgroup_border(colour = 'black', size = 2) +
            treemapify::geom_treemap_text(colour = 'black', place = 'topleft', reflow = TRUE, size = 16) +
            scale_fill_manual(values = palette, name = 'cluster') +
            labs(title = 'Clustered Gene Ontology Terms') + theme(legend.position = 'right')
        ggsave(s10_path, p_tree, width = 12, height = 8)
        engine <- 'treemapify'
    }
    rows[[length(rows) + 1]] <- manifest_row('FigureS10', 'treemap', s10_path, note = glue('{nrow(treemap_data)} ApoB GO terms by semantic cluster ({engine}); {cluster_note}'))
    say('Figures S9 and S10 done')
    do.call(rbind, rows)
}
