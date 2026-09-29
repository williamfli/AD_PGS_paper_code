# src/fig3_apoe.R -- Figure 3 (the APOE dissection of five PGS) and Figure S8 / Table S7.
#   make_figure3()   left column (a,c,e,g,i,k): effect [95% CI] of each PGS on six phenotypes,
#                    full PGS vs PGS without the APOE region; (k) the APOE gradient alone.
#                    right column (b,d,f,h,j): variant-level snpnet weights by genomic position,
#                    shape = consequence group.
#   make_figureS8()  APOE attenuation of the significant Stage-2 associations (heatmap and
#                    histogram) with the neuritic-plaque variance decomposition; Also Table S7.
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(grid)
    library(ComplexHeatmap)
    library(circlize)
})

# ---- Constants -----------------------------------------------------------------------
APOE_REGION_HG19 <- c(start = 45176340, end = 45447221)
APOE_CHROM       <- 19
FIG3_N_TOP_VARIANTS <- 5
FIG3_RIGHT_PANELS   <- c('b', 'd', 'f', 'h', 'j')

FIG3_SOURCE_COLOURS <- c('Full PGS' = '#E64B35', 'PGS without APOE region' = '#4DBBD5')
FIG3_CHROM_COLOURS  <- c('0' = '#E69F00', '1' = '#56B4E9')
FIG3_SHAPES <- c('FDR < 0.001' = 16, '0.001 <= FDR < 0.01' = 17, '0.01 <= FDR < 0.05' = 18, '0.05 <= FDR' = 22)
FIG3_CSQ_SHAPES <- c('Intronic' = 1, 'PAVs' = 2, 'UTR' = 3, 'PTVs' = 4, 'PCVs' = 5, 'Others' = 6, 'Unclassified' = 7)

fig3_wrap_label <- function(x, width = 22) paste(strwrap(x, width = width), collapse = '\n')

# ---- Left column: effect profiles -----------------------------------------------------
fig3_effect_table <- function(effects, errors, fdr, row, source) {
    f <- as.numeric(fdr[row, FIG3_PHENOTYPES])
    band <- rep('0.05 <= FDR', length(f))
    band[f < 0.05]  <- '0.01 <= FDR < 0.05'
    band[f < 0.01]  <- '0.001 <= FDR < 0.01'
    band[f < 0.001] <- 'FDR < 0.001'
    data.frame(variable = factor(FIG3_PHENOTYPES, levels = FIG3_PHENOTYPES),
               effects = as.numeric(effects[row, FIG3_PHENOTYPES]),
               errors = as.numeric(errors[row, FIG3_PHENOTYPES]),
               source = source,
               significance = factor(band, levels = names(FIG3_SHAPES)),
               alpha = factor(ifelse(f <= FDR_STAGE2, '1', '0.4'), levels = c('0.4', '1')),
               stringsAsFactors = FALSE)
}

# Scales use drop = FALSE and show.legend = TRUE so patchwork can collect one shared legend.
fig3_effect_panel <- function(df, ylab, colours, dodge_width, show_x = TRUE, colour_guide = TRUE) {
    df$source <- factor(df$source, levels = names(colours))
    p <- ggplot(df, aes(x = variable, y = effects, color = source, group = source, shape = significance, alpha = alpha)) +
        geom_point(size = 3, position = position_dodge(width = dodge_width), show.legend = TRUE) +
        geom_errorbar(aes(ymin = effects - 1.96 * errors, ymax = effects + 1.96 * errors),
                      width = 0.4, position = position_dodge(width = dodge_width), show.legend = TRUE) +
        geom_hline(yintercept = 0, linetype = 'dotted', color = 'black') +
        scale_alpha_manual(values = c('0.4' = 0.4, '1' = 1), drop = FALSE, labels = c('0.4' = 'FDR >= 0.1', '1' = 'FDR < 0.1')) +
        scale_shape_manual(values = FIG3_SHAPES, drop = FALSE) +
        scale_color_manual(values = colours, drop = FALSE, guide = if (colour_guide) 'legend' else 'none') +
        ylim(-0.4, 0.4) +
        theme_bw(base_size = 20) +
        theme(axis.text.y = element_text(angle = 0, hjust = 1),
              axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
              panel.grid.major.y = element_line(color = 'grey', linewidth = 0.5),
              panel.grid.minor.y = element_line(color = 'grey', linewidth = 0.25),
              plot.margin = unit(c(0, 0, 0, 0), 'cm'),
              plot.title = element_text(hjust = 0.5)) +
        labs(y = fig3_wrap_label(ylab), x = 'Phenotype (ROSMAP variable)', color = 'APOE Region Inclusion',
             shape = 'Significance', alpha = 'Opacity')
    if (!show_x) p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    p
}

fig3_effect_profiles <- function() {
    s2 <- load_stage2()
    na <- load_stage2_no_apoe()
    ao <- load_apoe_only()
    panels <- lapply(FIG3_PGS, function(pgs_name) {
        df <- rbind(fig3_effect_table(s2$effects, s2$errors, s2$fdr, pgs_name, 'Full PGS'),
                    fig3_effect_table(na$effects, na$errors, na$fdr, pgs_name, 'PGS without APOE region'))
        fig3_effect_panel(df, ylab = pgs_name, colours = FIG3_SOURCE_COLOURS, dodge_width = 0.4, show_x = FALSE)
    })
    apoe_df <- fig3_effect_table(ao$effects, ao$errors, ao$fdr, 'apoe_genotype', 'APOE region only')
    panels[[length(panels) + 1]] <- fig3_effect_panel(apoe_df, ylab = 'APOE region only', colours = c('APOE region only' = 'red'),
                                                      dodge_width = 0, show_x = TRUE, colour_guide = FALSE)
    wrap_plots(panels, ncol = 1) + plot_layout(guides = 'collect')
}

# ---- Right column: variant-level PGS profiles -----------------------------------------
pgs_model_file <- function(pgs_id) {
    dir <- optional_dir('pgs_models_dir')
    if (is.null(dir)) NA_character_ else glue('{dir}{pgs_id}.snpnetBETAs.tsv.gz')
}

load_chromosome_lengths <- function() memo('chrom_lengths', {
    tab <- read.delim(input_path('chromosome_lengths'), header = TRUE, sep = '\t', fill = TRUE)
    lengths <- as.numeric(gsub(',', '', tab$Length))
    chr_to_num <- function(chr) {
        chr <- as.character(chr)
        n <- suppressWarnings(as.numeric(chr))
        n[chr %in% c('X', 'XY')] <- 23
        n[chr == 'Y']  <- 24
        n[chr == 'MT'] <- 25
        n
    }
    offsets <- c(0, cumsum(lengths)[-length(lengths)])
    list(lengths = lengths, chr_to_num = chr_to_num, chr_to_pos = function(chr) offsets[chr_to_num(chr)], genome_end = sum(lengths))
})

load_variant_annotation <- function() memo('variant_annotation', {
    file <- optional_path('variant_annotation')
    say('Figure 3: reading variant annotation {file}')
    ann <- read.delim(file, header = TRUE, sep = '\t', fill = TRUE, quote = '', check.names = FALSE)
    rs_col <- if ('rsID' %in% names(ann)) 'rsID' else if ('ID_UKB' %in% names(ann)) 'ID_UKB' else
        stop('variant_annotation needs an rsID (or ID_UKB) column')
    if (!all(c('ID', 'Csq_group') %in% names(ann))) stop('variant_annotation needs ID and Csq_group columns')
    parts <- strsplit(as.character(ann$ID), ':')
    key <- vapply(parts, function(x) paste(ifelse(x[1] == 'XY', 'X', x[1]), x[2], sep = ':'), character(1))
    keep <- !duplicated(key, fromLast = TRUE)
    rs <- as.character(ann[[rs_col]][keep]); rs[is.na(rs)] <- ''
    list(rsid = setNames(rs, key[keep]), csq = setNames(as.character(ann$Csq_group[keep]), key[keep]))
})

read_pgs_model <- function(pgs_id) {
    chrom <- load_chromosome_lengths()
    ann <- load_variant_annotation()
    model <- read.delim(pgs_model_file(pgs_id), header = TRUE, sep = '\t', fill = TRUE, quote = '')
    id_col <- if ('X.ID' %in% names(model)) 'X.ID' else 'ID'
    parts <- strsplit(as.character(model[[id_col]]), ':')
    model$chromosome <- vapply(parts, `[`, character(1), 1)
    pos <- as.numeric(vapply(parts, `[`, character(1), 2))
    model$chrom_pos <- paste(ifelse(model$chromosome == 'XY', 'X', model$chromosome), pos, sep = ':')
    known <- model$chrom_pos %in% names(ann$csq)
    rs <- ann$rsid[model$chrom_pos]
    model$rsid <- ifelse(known & !is.na(rs) & nzchar(rs), rs, model$chrom_pos)
    model$csq_group <- ifelse(known, ann$csq[model$chrom_pos], 'Unclassified')
    model$csq_group[!(model$csq_group %in% names(FIG3_CSQ_SHAPES))] <- 'Unclassified'
    say('Figure 3: {pgs_id}: {nrow(model)} variants, {sum(!known)} not in the annotation, ',
        '{sum(model$csq_group == "Unclassified")} unclassified, {sum(!grepl("^rs", model$rsid))} without rsID')

    model$position <- pos + chrom$chr_to_pos(model$chromosome)
    model$chrom_parity <- factor(chrom$chr_to_num(model$chromosome) %% 2, levels = c('0', '1'))
    model$csq_group <- factor(model$csq_group, levels = names(FIG3_CSQ_SHAPES))
    model$top_snps <- FALSE
    model$top_snps[order(-abs(model$BETA))[seq_len(min(FIG3_N_TOP_VARIANTS, nrow(model)))]] <- TRUE
    model
}

fig3_pgs_profile_plot <- function(model, pgs_name) {
    chrom <- load_chromosome_lengths()
    apoe_x <- chrom$chr_to_pos(APOE_CHROM) + APOE_REGION_HG19
    ymax <- max(abs(model$BETA))
    labels_df <- model[model$top_snps, ]
    label_layer <- if (requireNamespace('ggrepel', quietly = TRUE)) {
        ggrepel::geom_text_repel(data = labels_df, aes(label = rsid), size = 6, show.legend = FALSE, min.segment.length = 0, seed = SPLIT_SEED)
    } else {
        geom_text(data = labels_df, aes(label = rsid), vjust = 1.5, size = 6, show.legend = FALSE)
    }
    ggplot(model, aes(x = position, y = BETA, color = chrom_parity, shape = csq_group)) +
        geom_point(size = 4) +
        scale_color_manual(values = FIG3_CHROM_COLOURS, guide = 'none') +
        scale_shape_manual(values = FIG3_CSQ_SHAPES, drop = TRUE) +
        theme_bw(base_size = 24) +
        labs(x = 'Genomic position (chromosome)', y = paste0(pgs_name, ' PGS weights'), shape = 'CSQ Group') +
        geom_hline(yintercept = 0, col = 'black') +
        geom_vline(xintercept = apoe_x, linetype = 'dashed', color = 'red', linewidth = 0.1) +
        scale_x_continuous(breaks = chrom$chr_to_pos(1:23) + chrom$lengths[1:23] / 2,
                           labels = c(1:9, ' ', 11, ' ', 13, ' ', 15, ' ', 17, ' ', 19, ' ', 21, ' ', 'X'),
                           limits = c(0, chrom$genome_end)) +
        theme(axis.text.x = element_text(angle = 90, size = 24), axis.text.y = element_text(size = 24),
              axis.title.x = element_text(size = 24), axis.title.y = element_text(size = 24),
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
        coord_cartesian(ylim = c(-ymax, ymax)) +
        label_layer
}

make_figure3 <- function() {
    rows <- list()

    say('Figure 3: effect profiles (a,c,e,g,i,k)')
    left_path <- fig_path('Figure3', 'Fig3_acegik_effect_profiles.pdf')
    ggsave(left_path, plot = fig3_effect_profiles(), width = 10, height = 5.5 * (length(FIG3_PGS) + 1), units = 'in', dpi = 300, limitsize = FALSE)
    rows[[1]] <- manifest_row('Figure3', 'a,c,e,g,i,k', left_path,
                              note = 'effect [95% CI] of full vs APOE-excluded PGS on six phenotypes; k = APOE gradient alone')

    # Right column: one panel per PGS when the model weights and annotation are provided.
    for (i in seq_along(FIG3_PGS)) {
        pgs_id <- FIG3_PGS_IDS[i]; pgs_name <- FIG3_PGS[i]; panel <- FIG3_RIGHT_PANELS[i]
        path <- fig_path('Figure3', glue('Fig3{panel}_{pgs_id}_variant_coefficients.pdf'))
        model_file <- pgs_model_file(pgs_id)
        missing <- c('pgs_models_dir', 'variant_annotation')[c(is.na(model_file) || !file.exists(model_file), is.null(optional_path('variant_annotation')))]
        if (length(missing) > 0) {
            rows[[length(rows) + 1]] <- manifest_row('Figure3', panel, path, 'skipped',
                                                     glue('PGS profile for {pgs_name} not built; set {paste(missing, collapse = " and ")} in scripts/paths.R'))
            next
        }
        say('Figure 3{panel}: PGS profile for {pgs_id} ({pgs_name})')
        model <- read_pgs_model(pgs_id)
        ggsave(path, plot = fig3_pgs_profile_plot(model, pgs_name), width = 25, height = 7)
        rows[[length(rows) + 1]] <- manifest_row('Figure3', panel, path,
                                                 note = glue('variant-level snpnet weights of {pgs_name} ({pgs_id}); {nrow(model)} variants, top {FIG3_N_TOP_VARIANTS} labelled'))
    }
    do.call(rbind, rows)
}

# ============================ Figure S8 + Table S7 ====================================
S8_FONT_SIZE <- 14
S8_ASTERISK_SIZE <- 9
COL_BLUE <- '#3C5488'
COL_RED <- '#E64B35'
COL_ATTENUATION_MID <- '#D9A3A2'
COL_ATTENUATION_HIGH <- '#B24745'
COL_NOT_TESTED <- 'grey75'
COL_VARIANCE_NON_APOE <- '#4DBBD5'
COL_VARIANCE_APOE <- '#3C5488'
WITH_APOE_SUFFIX <- '__with_apoe'

fdr_asterisks <- function(fdr) ifelse(is.na(fdr), '', ifelse(fdr < 0.001, '***', ifelse(fdr < 0.01, '**', ifelse(fdr < 0.05, '*', ''))))

# 'APOE' in italics inside a label; as.expression() so ComplexHeatmap's Legend() evaluates it.
italic_apoe <- function(before = '', after = '') as.expression(bquote(.(before) * italic('APOE') * .(after)))

s8_theme <- function() {
    theme_bw(base_size = S8_FONT_SIZE) +
        theme(text = element_text(size = S8_FONT_SIZE, colour = 'black'),
              plot.title = element_text(size = S8_FONT_SIZE, colour = 'black'),
              axis.title = element_text(size = S8_FONT_SIZE, colour = 'black'),
              axis.text = element_text(size = S8_FONT_SIZE, colour = 'black'),
              legend.title = element_text(size = S8_FONT_SIZE, colour = 'black'),
              legend.text = element_text(size = S8_FONT_SIZE, colour = 'black'))
}

# ---- Variance decomposition: R^2 each version of a PGS adds over the covariates ------------
build_decomposition_design <- function() {
    cohort <- load_cohort()
    excluded <- read_scores('pgs_scores_no_apoe')
    score_columns <- setdiff(colnames(excluded), '#projid')
    if (!all(score_columns %in% colnames(cohort))) stop('build_decomposition_design: APOE-excluded score columns not all present in the cohort table')
    idx <- match(cohort[[ID_COL]], excluded[['#projid']])
    if (anyNA(idx)) stop(glue('build_decomposition_design: {sum(is.na(idx))} cohort individuals missing from pgs_scores_no_apoe'))
    design <- cohort[, setdiff(colnames(cohort), c(ID_COL, score_columns)), drop = FALSE]
    design[score_columns] <- excluded[idx, score_columns]
    design[paste0(score_columns, WITH_APOE_SUFFIX)] <- cohort[, score_columns]
    for (column in colnames(design)) {
        if (is.numeric(design[[column]])) design[[column]][is.na(design[[column]])] <- mean(design[[column]], na.rm = TRUE)
    }
    design
}

r_squared <- function(formula, data) {
    fit <- lm(formula, data = data)
    response <- fit$model[[1]]
    1 - sum(residuals(fit)^2) / sum((response - mean(response))^2)
}

decompose_one_pair <- function(pgs, phenotype, design, covariate_r2_cache) {
    zscore <- function(x) (x - mean(x)) / sd(x)
    covariate_formula <- paste(COVARIATES, collapse = ' + ')
    model_data <- cbind(data.frame(y = zscore(design[[phenotype]]),
                                   pgs_excluded = zscore(design[[pgs]]),
                                   pgs_included = zscore(design[[paste0(pgs, WITH_APOE_SUFFIX)]])),
                        design[, COVARIATES])
    if (!exists(phenotype, envir = covariate_r2_cache))
        assign(phenotype, r_squared(as.formula(glue('y ~ {covariate_formula}')), model_data), envir = covariate_r2_cache)
    baseline_r2 <- get(phenotype, envir = covariate_r2_cache)
    var_excluded <- r_squared(as.formula(glue('y ~ pgs_excluded + {covariate_formula}')), model_data) - baseline_r2
    var_included <- r_squared(as.formula(glue('y ~ pgs_included + {covariate_formula}')), model_data) - baseline_r2
    data.frame(var_explained_pgs_excluded = var_excluded, var_explained_pgs_included = var_included,
               apoe_region_contribution = var_included - var_excluded)
}

decompose_pairs <- function(pgs, phenotype, design) {
    pairs <- data.frame(pgs = pgs, phenotype = phenotype, stringsAsFactors = FALSE)
    cache <- new.env(parent = emptyenv())
    decomposed <- do.call(rbind, Map(function(p, ph) decompose_one_pair(p, ph, design, cache), pairs$pgs, pairs$phenotype))
    result <- cbind(pairs, decomposed)
    stopifnot(isTRUE(all.equal(result$var_explained_pgs_included, result$var_explained_pgs_excluded + result$apoe_region_contribution)))
    result
}

# ---- Table S7: one row per association significant with the full PGS ---------------------
# Attenuation = 1 - clamp(beta_without / beta_with, 0, 1).
S7_EXPORT_COLUMNS <- c(
    pgs_name = 'UKB PGS',
    phenotype_name = 'Observed ROSMAP Phenotype',
    effect_with_apoe_region = 'Effect Size (with APOE region)',
    fdr_with_apoe_region = 'FDR (with APOE region)',
    effect_no_apoe_region = 'Effect Size (APOE region excluded)',
    fdr_no_apoe_region = 'FDR (APOE region excluded)',
    apoe_attenuation = 'APOE Effect Size Attenuation',
    var_explained_pgs_excluded = 'Variance Explained by APOE-excluded PGS',
    var_explained_pgs_included = 'Variance Explained by APOE-included PGS',
    apoe_region_contribution = 'APOE-region Variance Contribution')

apoe_variance_decomposition_table <- function(design) {
    with_apoe <- lapply(load_stage2(mapped = FALSE)[c('effects', 'fdr')], as.matrix)
    without_apoe <- lapply(load_stage2_no_apoe(mapped = FALSE)[c('effects', 'fdr')], function(m)
        as.matrix(m)[rownames(with_apoe$effects), colnames(with_apoe$effects)])

    significant_cells <- which(with_apoe$fdr < FDR_STAGE2, arr.ind = TRUE)
    pairs <- data.frame(pgs = rownames(with_apoe$effects)[significant_cells[, 'row']],
                        phenotype = colnames(with_apoe$effects)[significant_cells[, 'col']], stringsAsFactors = FALSE)
    cells <- function(m) m[cbind(pairs$pgs, pairs$phenotype)]
    pairs$effect_with_apoe_region <- cells(with_apoe$effects)
    pairs$fdr_with_apoe_region <- cells(with_apoe$fdr)
    pairs$effect_no_apoe_region <- cells(without_apoe$effects)
    pairs$fdr_no_apoe_region <- cells(without_apoe$fdr)
    pairs$effect_attenuation_ratio <- pairs$effect_no_apoe_region / pairs$effect_with_apoe_region
    pairs$apoe_attenuation <- 1 - pmin(pmax(pairs$effect_attenuation_ratio, 0), 1)

    decomposition <- decompose_pairs(pairs$pgs, pairs$phenotype, design)
    pairs <- cbind(pairs, decomposition[, c('var_explained_pgs_excluded', 'var_explained_pgs_included', 'apoe_region_contribution')])
    stopifnot(!any(is.na(pairs[, sapply(pairs, is.numeric)])))
    pairs$pgs_name <- mapper_pgs(pairs$pgs)
    pairs$phenotype_name <- mapper_pheno(pairs$phenotype)

    output_table <- pairs[, names(S7_EXPORT_COLUMNS)]
    colnames(output_table) <- unname(S7_EXPORT_COLUMNS)
    output_table
}

make_figureS8 <- function() {
    # ---- Table S7 --------------------------------------------------------------------
    design <- build_decomposition_design()
    table_s7 <- apoe_variance_decomposition_table(design)
    s7_path <- table_path('TableS7_apoe_variance_decomposition.tsv')
    write.table(table_s7, file = s7_path, sep = '\t', row.names = FALSE, quote = FALSE)
    say('Table S7: {nrow(table_s7)} significant associations, ',
        '{sum(table_s7[["FDR (APOE region excluded)"]] < FDR_STAGE2)} still significant without APOE')

    # ---- Attenuation matrix (phenotypes x PGS, as in Figure 2) -------------------------
    s2 <- load_stage2(mapped = FALSE)
    s2na <- load_stage2_no_apoe(mapped = FALSE)
    effects_with <- as.matrix(s2$effects)
    fdr_with <- as.matrix(s2$fdr)[rownames(effects_with), colnames(effects_with)]
    effects_without <- as.matrix(s2na$effects)[rownames(effects_with), colnames(effects_with)]

    retained_ratio <- t(effects_without / effects_with)
    retained_ratio[t(fdr_with) >= FDR_STAGE2] <- NA
    attenuation <- 1 - pmin(pmax(retained_ratio, 0), 1)
    asterisks <- t(fdr_asterisks(fdr_with))
    rownames(attenuation) <- mapper_pheno(rownames(attenuation))
    colnames(attenuation) <- mapper_pgs(colnames(attenuation))
    dimnames(asterisks) <- dimnames(attenuation)
    dendrograms <- stage2_dendrograms()

    # ---- Variance decomposition of neuritic plaque burden for all 12 PGS ----------------
    example_phenotype_id <- colnames(effects_with)[mapper_pheno(colnames(effects_with)) == 'Neuritic plaque burden (5 regions)']
    stopifnot(length(example_phenotype_id) == 1)
    example <- decompose_pairs(rownames(effects_with), example_phenotype_id, design)
    example$pgs_name <- mapper_pgs(example$pgs)
    annotation_values <- as.matrix(example[match(colnames(attenuation), example$pgs_name),
                                           c('apoe_region_contribution', 'var_explained_pgs_excluded')])
    rownames(annotation_values) <- colnames(attenuation)
    stopifnot(!any(is.na(annotation_values)))

    fs <- S8_FONT_SIZE
    variance_annotation <- HeatmapAnnotation(
        plaque_variance = anno_barplot(annotation_values, gp = gpar(fill = c(COL_VARIANCE_APOE, COL_VARIANCE_NON_APOE), col = NA),
                                       bar_width = 0.75, height = unit(5.4, 'cm'), axis_param = list(gp = gpar(fontsize = fs))),
        which = 'column', annotation_name_side = 'left',
        annotation_label = list(plaque_variance = expression(paste('Neuritic plaque variance (', R^2, ')'))),
        annotation_name_gp = gpar(fontsize = fs, col = 'black'), annotation_name_rot = 0)

    joint_heatmap <- Heatmap(
        attenuation, name = 'APOE attenuation',
        col = colorRamp2(c(0, 0.5, 1), c('white', COL_ATTENUATION_MID, COL_ATTENUATION_HIGH)),
        na_col = COL_NOT_TESTED, rect_gp = gpar(col = 'grey85', lwd = 0.4),
        cluster_rows = dendrograms$row, cluster_columns = dendrograms$col,
        top_annotation = variance_annotation,
        row_dend_side = 'right', row_dend_width = unit(1.2, 'cm'), column_dend_height = unit(1.2, 'cm'),
        row_names_side = 'left', column_names_side = 'bottom', column_names_rot = 90,
        row_names_max_width = max_text_width(rownames(attenuation), gp = gpar(fontsize = fs)),
        column_names_max_height = max_text_width(colnames(attenuation), gp = gpar(fontsize = fs)),
        row_title = 'Observed ROSMAP Phenotypes', column_title = 'UKB PGS',
        row_title_side = 'left', column_title_side = 'bottom',
        row_title_gp = gpar(fontsize = fs, col = 'black'), column_title_gp = gpar(fontsize = fs, col = 'black'),
        width = ncol(attenuation) * unit(7, 'mm'), height = nrow(attenuation) * unit(7, 'mm'),
        row_names_gp = gpar(fontsize = fs, col = 'black'), column_names_gp = gpar(fontsize = fs, col = 'black'),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(asterisks[i, j], x, y, gp = gpar(fontsize = S8_ASTERISK_SIZE, col = 'black'))
        },
        heatmap_legend_param = list(title = italic_apoe(after = ' attenuation'), legend_direction = 'vertical',
                                    at = c(0, 0.25, 0.5, 0.75, 1),
                                    title_gp = gpar(fontsize = fs, col = 'black'), labels_gp = gpar(fontsize = fs, col = 'black')))

    variance_legend <- Legend(labels = c(italic_apoe(after = '-excluded PGS'), italic_apoe(after = ' region')),
                              legend_gp = gpar(fill = c(COL_VARIANCE_NON_APOE, COL_VARIANCE_APOE)),
                              title = 'Variance explained in neuritic plaque burden', direction = 'horizontal',
                              title_gp = gpar(fontsize = fs, col = 'black'), labels_gp = gpar(fontsize = fs, col = 'black'))

    s8a_path <- fig_path('FigureS8', 'FigS8a_apoe_attenuation_with_variance_decomposition.pdf')
    pdf(s8a_path, width = 11.1, height = 19.4)
    draw(joint_heatmap, heatmap_legend_side = 'right', annotation_legend_list = list(variance_legend),
         annotation_legend_side = 'top', padding = unit(c(2, 2, 2, 2), 'mm'))
    invisible(dev.off())

    # ---- S8b: histogram of the attenuation values ----------------------------------------
    attenuation_values <- attenuation[!is.na(attenuation)]
    histogram_data <- data.frame(
        attenuation = pmin(attenuation_values, 1 - 1e-9),
        kind = factor(ifelse(attenuation_values == 0, 'none', 'attenuated'), levels = c('none', 'attenuated')))
    kind_counts <- table(histogram_data$kind)
    histogram <- ggplot(histogram_data, aes(attenuation, fill = kind)) +
        geom_histogram(breaks = seq(0, 1, by = 0.05), closed = 'left', colour = 'grey30', linewidth = 0.2) +
        scale_x_continuous(limits = c(-0.02, 1.02), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        scale_fill_manual(values = c(none = COL_BLUE, attenuated = COL_RED),
                          labels = c(italic_apoe('No ', glue('-locus variant (n={kind_counts[["none"]]})')),
                                     as.expression(glue('Attenuated (n={kind_counts[["attenuated"]]})')))) +
        labs(x = italic_apoe(after = ' attenuation'), y = expression(PGS %*% phenotype ~ associations), fill = NULL) +
        s8_theme() + theme(legend.position = 'top', panel.grid.minor = element_blank())
    s8b_path <- fig_path('FigureS8', 'FigS8b_apoe_attenuation_histogram.pdf')
    ggsave(s8b_path, histogram, width = 6, height = 5.5, dpi = 300)

    rbind(
        manifest_row('FigureS8', 'a', s8a_path, note = glue('APOE attenuation heatmap in Figure 2 order with neuritic-plaque variance bars; {sum(!is.na(attenuation))} significant cells')),
        manifest_row('FigureS8', 'b', s8b_path, note = glue('attenuation histogram: {kind_counts[["none"]]} no APOE-locus variant, {kind_counts[["attenuated"]]} attenuated')),
        manifest_row('TableS7', 'table', s7_path, note = glue('{nrow(table_s7)} rows')))
}
