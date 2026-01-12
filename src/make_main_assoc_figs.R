library(ComplexHeatmap)
library(circlize)
library(grid)
library(glue)
create_scatter_plot <- function(merged_data, pgs_name, pheno_name, 
                                 effects, errors, fdr, mapper_pgs, mapper_pheno, filename) {
    pdf(filename, width = 8, height = 16)
    par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
    all_x_vals <- vector("list", 2)
    all_y_vals <- vector("list", 2)
    for (i in 1:2) {
        for (j in 1:2) {
            x_vals <- merged_data[[pgs_name[i]]]
            y_vals <- merged_data[[pheno_name[j]]]
            valid_idx <- !is.na(x_vals) & !is.na(y_vals)
            x_vals <- x_vals[valid_idx]
            y_vals <- y_vals[valid_idx]
            x_vals <- (x_vals - mean(x_vals)) / sd(x_vals)
            all_x_vals[[i]] <- c(all_x_vals[[i]], x_vals)
            all_y_vals[[j]] <- c(all_y_vals[[j]], y_vals)
        }
    }
    x_lim <- lapply(all_x_vals, range)
    y_lim <- lapply(all_y_vals, range)
    for (j in 1:2) {
        for (i in 1:2) {
            x_vals <- merged_data[[pgs_name[i]]]
            y_vals <- merged_data[[pheno_name[j]]]
            valid_idx <- !is.na(x_vals) & !is.na(y_vals)
            x_vals <- x_vals[valid_idx]
            y_vals <- y_vals[valid_idx]
            x_vals <- (x_vals - mean(x_vals)) / sd(x_vals)
            mapped_pgs <- mapper_pgs(pgs_name[i])
            mapped_pheno <- mapper_pheno(pheno_name[j])
            slope <- as.numeric(effects[mapped_pgs, mapped_pheno])
            slope_plot <- slope * sd(y_vals)
            plot(x_vals, y_vals, 
                 xlim = x_lim[[i]], ylim = y_lim[[j]],
                 xlab = mapped_pgs, 
                 ylab = mapped_pheno,
                 main = glue("{mapped_pgs} vs {mapped_pheno}"),
                 pch = 19, col = rgb(0, 0, 1, 0.3), cex = 1.,
                 cex.lab = 1.8, cex.main = 1.5, cex.axis = 1.5,)
            se <- as.numeric(errors[mapped_pgs, mapped_pheno])
            upper_slope <- (slope + 1.96 * se) * sd(y_vals)
            lower_slope <- (slope - 1.96 * se) * sd(y_vals)
            x_range <- x_lim[[i]]
            y_line <- slope_plot * x_range + mean(y_vals)
            y_upper <- upper_slope * x_range + mean(y_vals)
            y_lower <- lower_slope * x_range + mean(y_vals)
            polygon(c(x_range, rev(x_range)), c(y_upper, rev(y_lower)),
                    col = rgb(1, 0, 0, 0.2), border = NA)
            lines(x_range, y_line, col = "red", lwd = 3)
            fdr_val <- as.numeric(fdr[mapped_pgs, mapped_pheno])
            legend("topright", 
                   legend = c(glue("Beta: {round(slope, 3)} ± {round(se, 3)}"),
                             glue("FDR: {format(fdr_val, scientific = TRUE, digits = 2)}")),
                   col = "red", lwd = 3, cex = 1.2, bty = "n")
        }
    }
    dev.off()

}


create_heatmap <- function(effects, fdr, filename) {
    heatmap_matrix <- as.matrix(effects)
    asterisks_matrix <- matrix("", nrow = nrow(fdr), ncol = ncol(fdr))
    asterisks_matrix[fdr < 0.05 & fdr >= 0.01] <- "*"
    asterisks_matrix[fdr < 0.01 & fdr >= 0.001] <- "**"
    asterisks_matrix[fdr < 0.001] <- "***"
    # asterisks_matrix[fdr < 0.5] <- '*'
    rownames(asterisks_matrix) <- rownames(fdr)
    colnames(asterisks_matrix) <- colnames(fdr)
    set.seed(123)
    heatmap_matrix <- t(heatmap_matrix)
    asterisks_matrix <- t(asterisks_matrix)
    p <- Heatmap(heatmap_matrix,
                 name = "Effect size",
                 width = ncol(heatmap_matrix)*unit(12, "mm"), 
                 height = nrow(heatmap_matrix)*unit(12, "mm"),
                 column_dend_height = unit(2, "cm"), 
                 row_dend_width = unit(2, "cm"),
                 col = colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "orange")),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(asterisks_matrix[i, j], x, y, gp = gpar(fontsize = 10))
                 },
                 row_title_gp = gpar(fontsize = 26),
                 cluster_rows = TRUE,
                 cluster_columns = TRUE,
                 show_row_names = TRUE,
                 show_column_names = TRUE,
                 row_names_gp = gpar(fontsize = 26),
                 column_names_gp = gpar(fontsize = 26),
                 heatmap_legend_param = list(
                     title = "Effect size (beta)", 
                     legend_direction = "horizontal",
                     title_gp = gpar(fontsize = 26),
                     labels_gp = gpar(fontsize = 26)
                 )
    )

    # Save row and column dendrograms as R objects
    row_dend <- row_dend(p)
    col_dend <- column_dend(p)
    saveRDS(row_dend, file = glue("{filename}_row_dend.rds"))
    saveRDS(col_dend, file = glue("{filename}_col_dend.rds"))

    # Save as PDF
    pdf(glue('{filename}.pdf'), width=100, height = 50)
    draw(p,  heatmap_legend_side = "top",  padding = unit(c(4, 4, 4, 4), "cm"))
    dev.off()
}

create_individuals_heatmap <- function(phenotypes, pgs, mapper_pheno, mapper_pgs, full_input, pgs_vars, train_idx, file) {
    rownames(phenotypes) <- phenotypes[["X.projid"]]
    phenotypes <- phenotypes[, colnames(phenotypes) != "X.projid"]
    # Flip the "Consensus cognitive diagnosis" variable
    if ("cogdx" %in% colnames(phenotypes)) {
        phenotypes$cogdx <- -phenotypes$cogdx
    }
    heatmap_matrix <- as.matrix(phenotypes)
    colnames(heatmap_matrix) <- sapply(colnames(heatmap_matrix), mapper_pheno)
    # Z-transform (standardize) each column, removing NAs
    z_transform <- function(x) {
        mu <- mean(x, na.rm = TRUE)
        sigma <- sd(x, na.rm = TRUE)
        (x - mu) / sigma
    }
    heatmap_matrix <- apply(heatmap_matrix, 2, z_transform)
    cat("Max of heatmap_matrix:", max(heatmap_matrix, na.rm = TRUE), "\n")
    cat("Min of heatmap_matrix:", min(heatmap_matrix, na.rm = TRUE), "\n")


    # Impute missing values with column means for clustering
    impute_col_means <- function(x) {
        inds <- is.na(x)
        x[inds] <- mean(x, na.rm = TRUE)
        x
    }
    heatmap_matrix_imputed <- apply(heatmap_matrix, 2, impute_col_means)

    heatmap_matrix <- t(heatmap_matrix)
    heatmap_matrix_imputed <- t(heatmap_matrix_imputed)

    p <- Heatmap(heatmap_matrix,
                width = ncol(heatmap_matrix)*unit(.24, "mm"), 
                height = nrow(heatmap_matrix)*unit(12, "mm"),
                column_dend_height = unit(2, "cm"), 
                row_dend_width = unit(2, "cm"),
                col = colorRamp2(c(-2, 0, 2), c("violet", "white", "red")),
                row_title = "ROSMAP Phenotypes",
                column_title = "Individuals",
                row_title_gp = gpar(fontsize = 22),
                column_title_gp = gpar(fontsize = 22),
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                clustering_distance_rows = dist(heatmap_matrix_imputed),
                clustering_distance_columns = dist(t(heatmap_matrix_imputed)),
                show_row_names = TRUE,
                show_column_names = FALSE,
                row_names_gp = gpar(fontsize = 22),
                column_names_gp = gpar(fontsize = 22),
                heatmap_legend_param = list(
                    title = "Z-score",
                    legend_direction = "horizontal",
                    title_gp = gpar(fontsize = 22),
                    labels_gp = gpar(fontsize = 22)
                )
    )


    pdf(file, width=50, height = 50)
    draw(p, heatmap_legend_side = "top", padding = unit(c(4, 4, 4, 4), "cm"))
    dev.off()


    # Prepare PGS matrix for plotting
    rownames(pgs) <- pgs[["X.projid"]]
    pgs <- pgs[, colnames(pgs) != "X.projid"]
    pgs_matrix <- as.matrix(pgs)
    colnames(pgs_matrix) <- sapply(colnames(pgs_matrix), mapper_pgs)

    # Z-transform (standardize) each column, removing NAs
    pgs_matrix <- apply(pgs_matrix, 2, z_transform)

    # Transpose to match heatmap orientation (rows = PGS, columns = individuals)
    pgs_matrix <- t(pgs_matrix)

    # Use the same column order as the phenotypes heatmap
    pgs_matrix <- pgs_matrix[, colnames(heatmap_matrix), drop = FALSE]

    # Plot PGS heatmap with a different color scheme
    # Only label selected PGS rows (after mapping)
    label_pgs_raw <- c(
        "cancer1044_SCORE1_AVG", "HC1584_SCORE1_AVG", "FH1044_SCORE1_AVG", "FH1263_SCORE1_AVG",
        "BIN_FC30006164_SCORE1_AVG", "BIN_FC40002267_SCORE1_AVG", "BIN22126_SCORE1_AVG",
        "INI30690_SCORE1_AVG", "INI30640_SCORE1_AVG", "INI30710_SCORE1_AVG",
        "INI30780_SCORE1_AVG", "INI25869_SCORE1_AVG"
    )
    # Subset pgs_matrix to only include rows in label_pgs after mapping
    label_pgs <- sapply(label_pgs_raw, mapper_pgs)
    pgs_matrix <- pgs_matrix[rownames(pgs_matrix) %in% label_pgs, , drop = FALSE]
    show_row_names_vec <- rownames(pgs_matrix) %in% label_pgs

    p_pgs <- Heatmap(
        pgs_matrix,
        width = ncol(pgs_matrix) * unit(.24, "mm"),
        height = nrow(pgs_matrix) * unit(7, "mm"),
        column_dend_height = unit(2, "cm"),
        row_dend_width = unit(2, "cm"),
        col = colorRamp2(c(-2, 0, 2), c("navy", "white", "gold")),
        row_title = "Polygenic Scores",
        column_title = "Individuals",
        row_title_gp = gpar(fontsize = 22),
        column_title_gp = gpar(fontsize = 22),
        cluster_rows = TRUE,
        cluster_columns = TRUE, # keep same order as phenotypes
        show_row_names = TRUE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 22),
        column_names_gp = gpar(fontsize = 22),
        heatmap_legend_param = list(
            title = "Z-score",
            legend_direction = "horizontal",
            title_gp = gpar(fontsize = 22),
            labels_gp = gpar(fontsize = 22)
        ),
        row_labels = ifelse(show_row_names_vec, rownames(pgs_matrix), "")
    )

    pdf(gsub(".pdf$", "_pgs.pdf", file), width = 35, height = 25)
    draw(p_pgs, heatmap_legend_side = "top", padding = unit(c(4, 4, 4, 4), "cm"))
    dev.off()

    create_individuals_heatmap_pca(full_input, pgs_vars, train_idx, gsub(".pdf$", "_pca.pdf", file), colnames(heatmap_matrix))
}

create_individuals_heatmap_pca <- function(full_input, pgs_vars, train_idx, file, colorder) {
  train_data <- full_input[train_idx, , drop = FALSE]
  pgs_train <- train_data[, pgs_vars]
  pca_result <- prcomp(pgs_train, scale. = TRUE)
  pca_all <- predict(pca_result, newdata = full_input[, pgs_vars, drop = FALSE])

  pca_matrix <- as.matrix(pca_all)
  colnames(pca_matrix) <- paste0("PC", 1:ncol(pca_matrix))
  rownames(pca_matrix) <- full_input$X.projid
  pca_matrix <- t(pca_matrix)
  pca_matrix <- pca_matrix[, colorder, drop = FALSE]

 p <- Heatmap(
    pca_matrix,
    width = ncol(pca_matrix) * unit(.24, "mm"),
    height = nrow(pca_matrix) * unit(7, "mm"),
    column_dend_height = unit(2, "cm"),
    row_dend_width = unit(2, "cm"),
    
    col = colorRamp2(c(-2, 0, 2), c("#0072B2", "white", "#D55E00")), # blue-white-orange (colorblind friendly)
    row_title = "PCA of Polygenic Scores",
    column_title = "Individuals",
    row_title_gp = gpar(fontsize = 22),
    column_title_gp = gpar(fontsize = 22),
    cluster_rows = TRUE,
    cluster_columns = TRUE, # keep same order as phenotypes
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 22),
    column_names_gp = gpar(fontsize = 22),
    heatmap_legend_param = list(
        title = "PGS PCs",
        legend_direction = "horizontal",
        title_gp = gpar(fontsize = 22),
        labels_gp = gpar(fontsize = 22)
    ))
    pdf(file, width = 35, height = 25)
    draw(p, heatmap_legend_side = "top", padding = unit(c(4, 4, 4, 4), "cm"))
    dev.off()
}