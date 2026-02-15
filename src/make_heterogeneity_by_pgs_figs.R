library(ggplot2)
library(scales)
library(glue)
library(ggfortify)
library(ggrepel)
library(xgboost)
library(Rtsne)
library(ComplexHeatmap)
library(circlize)
library(SHAPforxgboost)
set.seed(123)


pearson_correlation <- function(x, y) {
    return(cor(x, y, method = "pearson"))
}
kendall_correlation <- function(x, y) {
    return(cor(x, y, method = "kendall"))
}
compare_xgboost_models <- function(full_input, pgs_variables, selected_phenotype_variables, mapper_pgs, mapper_pheno, plot_directory) {
    covariates <- c("PC1_AVG","PC2_AVG","PC3_AVG","PC4_AVG","PC5_AVG","PC6_AVG","PC7_AVG","PC8_AVG","PC9_AVG","PC10_AVG",
                    "msex","age_death","age_bl","kronos")
    # Define scheme names for performance metrics:
    # Scheme 1: PCs (using PC1 through PC8 together)
    # Scheme 2: One PGS variable at a time (using columns in pgs_variables)
    # Scheme 3: PCs independently (each of PC1, PC2, …, PC8)
    # Scheme 4: All PGSs together (using the full set of pgs_variables)
    # Scheme 5: Only covariates
    scheme_names <- c("PCs", paste0("PC", 1:8), sapply(pgs_variables, mapper_pgs), "all_PGSs", "covariates")
    models <- data.frame(matrix(NA, nrow = length(selected_phenotype_variables), ncol = length(scheme_names)),
                         stringsAsFactors = FALSE)
    rownames(models) <- selected_phenotype_variables
    colnames(models) <- scheme_names
    
    params <- list(objective = "reg:squarederror", eval_metric = "rmse")

    # Preprocess once: impute missing values and perform train-test split
    # Split the data into training and testing sets before imputation
    n <- nrow(full_input)
    train_idx <- sample(seq_len(n), size = floor(0.7 * n))
    test_idx <- setdiff(seq_len(n), train_idx)
    
    train_data <- full_input[train_idx, , drop = FALSE]
    test_data  <- full_input[test_idx, , drop = FALSE]
    
    # Impute missing values separately for training data
    for (col in colnames(train_data)) {
        train_data[[col]][is.na(train_data[[col]])] <- mean(train_data[[col]], na.rm = TRUE)
    }
    
    # Impute missing values separately for testing data
    for (col in colnames(test_data)) {
        test_data[[col]][is.na(test_data[[col]])] <- mean(test_data[[col]], na.rm = TRUE)
    }
    
    # Recompute PCA on the imputed train PGS variables and extract PC1–PC8
    pgs_train <- train_data[, pgs_variables]
    pca_result <- prcomp(pgs_train, scale. = TRUE)
    pc_train <- as.matrix(pca_result$x[, 1:8])
    
    # Apply PCA transformation on imputed test PGS variables using the model from training data
    pgs_test <- test_data[, pgs_variables]
    pc_test <- predict(pca_result, newdata = pgs_test)[, 1:8]
    
    # Loop over each phenotype
    for (pheno in selected_phenotype_variables) {

        # Extract response variables from imputed train and test data
        y_train <- train_data[[pheno]]
        y_test  <- test_data[[pheno]]
        
        
        # ----- Scheme 1: Use PC1 through PC8 together plus covariates -----
        X_train <- cbind(pc_train, as.matrix(train_data[, covariates]))
        dtrain  <- xgb.DMatrix(data = X_train, label = y_train)
        model   <- xgb.train(params = params, data = dtrain, nrounds = 200)
      
        X_test  <- cbind(pc_test, as.matrix(test_data[, covariates]))
        preds   <- predict(model, newdata = X_test)
        r_value <- cor(preds, y_test)
        models[pheno, "PCs"] <- r_value
        # Save Shapley values for Scheme 1 (PCs model)
        shap_values <- shap.values(xgb_model = model, X_train = X_train)
        shap_df <- as.data.frame(shap_values$shap_score)
        colnames(shap_df) <- c(paste0("PC", 1:8), covariates)
        shap_means <- data.frame(t(colMeans(abs(shap_df[, paste0("PC", 1:8)]))))
        colnames(shap_means) <- paste0("PC", 1:8)
        write.table(shap_means,
            file = file.path(plot_directory, paste0("shapley_scores_PCs_", pheno, ".tsv")),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
        
        # ----- Scheme 2: Train xgboost on one PGS at a time plus covariates -----
        for (pg in pgs_variables) {
            X_train_pg <- as.matrix(cbind(train_data[, pg, drop = FALSE], train_data[, covariates]))
            dtrain_pg  <- xgb.DMatrix(data = X_train_pg, label = y_train)
            model_pg   <- xgb.train(params = params, data = dtrain_pg, nrounds = 200)
            
            X_test_pg <- as.matrix(cbind(test_data[, pg, drop = FALSE], test_data[, covariates]))
            preds_pg  <- predict(model_pg, newdata = X_test_pg)
            r_value_pg <- cor(preds_pg, y_test)
            models[pheno, mapper_pgs(pg)] <- r_value_pg
        }
        
        # ----- Scheme 3: Train xgboost using each PC independently plus covariates -----
        for (pc_index in 1:ncol(pc_train)) {
            pc_name <- paste0("PC", pc_index)
            X_train_pc <- as.matrix(cbind(pc_train[, pc_index, drop = FALSE], train_data[, covariates]))
            dtrain_pc  <- xgb.DMatrix(data = X_train_pc, label = y_train)
            model_pc   <- xgb.train(params = params, data = dtrain_pc, nrounds = 200)
            X_test_pc  <- as.matrix(cbind(pc_test[, pc_index, drop = FALSE], test_data[, covariates]))
            preds_pc   <- predict(model_pc, newdata = X_test_pc)
            r_value_pc <- cor(preds_pc, y_test)
            models[pheno, pc_name] <- r_value_pc
        }
        
        # ----- Scheme 4: Train xgboost on all PGSs together plus covariates -----
        X_train_all <- as.matrix(cbind(train_data[, pgs_variables], train_data[, covariates]))
        dtrain_all  <- xgb.DMatrix(data = X_train_all, label = y_train)
        model_all   <- xgb.train(params = params, data = dtrain_all, nrounds = 200)
        
        X_test_all  <- as.matrix(cbind(test_data[, pgs_variables], test_data[, covariates]))
        preds_all   <- predict(model_all, newdata = X_test_all)
        r_value_all <- cor(preds_all, y_test)
        models[pheno, "all_PGSs"] <- r_value_all
        # Save Shapley values for Scheme 4 (all PGSs model)
        shap_values_all <- shap.values(xgb_model = model_all, X_train = X_train_all)
        shap_df_all <- as.data.frame(shap_values_all$shap_score)
        colnames(shap_df_all) <- c(pgs_variables, covariates)
        shap_means_all <- data.frame(t(colMeans(abs(shap_df_all[, pgs_variables]))))
        colnames(shap_means_all) <- pgs_variables
        write.table(shap_means_all,
            file = file.path(plot_directory, paste0("shapley_scores_all_PGSs_", pheno, ".tsv")),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
        # ----- Scheme 5: Train xgboost using only the covariates -----
        X_train_cov <- as.matrix(train_data[, covariates])
        dtrain_cov  <- xgb.DMatrix(data = X_train_cov, label = y_train)
        model_cov   <- xgb.train(params = params, data = dtrain_cov, nrounds = 200)
        X_test_cov  <- as.matrix(test_data[, covariates])
        preds_cov   <- predict(model_cov, newdata = X_test_cov)
        r_value_cov <- cor(preds_cov, y_test)
        models[pheno, "covariates"] <- r_value_cov
        # ----- Scheme 6: Train xgboost using covariates + apoe_genotype -----
        covariates_apoe <- c(covariates, "apoe_genotype")
        X_train_apoe <- as.matrix(train_data[, covariates_apoe])
        dtrain_apoe  <- xgb.DMatrix(data = X_train_apoe, label = y_train)
        model_apoe   <- xgb.train(params = params, data = dtrain_apoe, nrounds = 200)
        X_test_apoe  <- as.matrix(test_data[, covariates_apoe])
        preds_apoe   <- predict(model_apoe, newdata = X_test_apoe)
        r_value_apoe <- cor(preds_apoe, y_test)
        models[pheno, "covariates_and_APOE"] <- r_value_apoe
    }

    # Create heatmap of mean Shapley contributions for all_PGSs model
    # Read in all Shapley score files and compute means
    shapley_means <- matrix(NA, nrow = length(selected_phenotype_variables), ncol = length(pgs_variables))
    rownames(shapley_means) <- selected_phenotype_variables
    colnames(shapley_means) <- pgs_variables
    for (pheno in selected_phenotype_variables) {
        shap_file <- file.path(plot_directory, paste0("shapley_scores_all_PGSs_", pheno, ".tsv"))
        if (file.exists(shap_file)) {
            shap_data <- read.delim(shap_file, header = TRUE, sep = "\t")
            pgs_shap <- shap_data[, pgs_variables, drop = FALSE]
            shapley_means[pheno, ] <- as.numeric(pgs_shap)
        }
    }
    # Normalize so that each row sums to 1
    shapley_means <- shapley_means / rowSums(shapley_means, na.rm = TRUE)
    rownames(shapley_means) <- sapply(rownames(shapley_means), mapper_pheno)
    colnames(shapley_means) <- sapply(colnames(shapley_means), mapper_pgs)
    ht_shapley <- Heatmap(
        shapley_means,
        name = "Normalized Mean |SHAP|",
        col = colorRamp2(c(0, max(shapley_means, na.rm = TRUE)), c("white", "darkred")),
        row_names_side = "left",
        column_names_side = "bottom",
        column_names_rot = 90,
        column_title = "PGS Variables",
        row_title = "Phenotypes",
        heatmap_legend_param = list(title = "Normalized Mean |SHAP|"),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_names_max_width = unit(8, "cm")
    )
    
    pdf(file.path(plot_directory, "shapley_mean_contributions_heatmap.pdf"), width = 12, height = 10)
    draw(ht_shapley)
    dev.off()
    
    # Create heatmap of mean Shapley contributions for PCs model
        shapley_means_pcs <- matrix(NA, nrow = length(selected_phenotype_variables), ncol = 8)
        rownames(shapley_means_pcs) <- selected_phenotype_variables
        colnames(shapley_means_pcs) <- paste0("PC", 1:8)
        for (pheno in selected_phenotype_variables) {
            shap_file <- file.path(plot_directory, paste0("shapley_scores_PCs_", pheno, ".tsv"))
            if (file.exists(shap_file)) {
                shap_data <- read.delim(shap_file, header = TRUE, sep = "\t")
                pcs_shap <- shap_data[, paste0("PC", 1:8), drop = FALSE]
                shapley_means_pcs[pheno, ] <- as.numeric(pcs_shap)
            }
        }
        # Normalize so that each row sums to 1
        shapley_means_pcs <- shapley_means_pcs / rowSums(shapley_means_pcs, na.rm = TRUE)
        rownames(shapley_means_pcs) <- sapply(rownames(shapley_means_pcs), mapper_pheno)
        ht_shapley_pcs <- Heatmap(
            shapley_means_pcs,
            name = "Normalized Mean |SHAP|",
            col = colorRamp2(c(0, max(shapley_means_pcs, na.rm = TRUE)), c("white", "darkorange")),
            row_names_side = "left",
            column_names_side = "bottom",
            column_names_rot = 90,
            column_title = "PGS Principal Components",
            row_title = "Phenotypes",
            heatmap_legend_param = list(title = "Normalized Mean |SHAP|"),
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            row_names_max_width = unit(8, "cm")
        )
        pdf(file.path(plot_directory, "shapley_mean_contributions_PCs_heatmap.pdf"), width = 12, height = 10)
        draw(ht_shapley_pcs)
        dev.off()
    rownames(models) <- sapply(rownames(models), mapper_pheno)
    write.table(models,
                file = file.path(plot_directory, "models.tsv"),
                sep = "\t",
                row.names = TRUE,
                quote = FALSE)
    return(list(
        train_idx = train_idx,
        test_idx = test_idx,
        pca_model = pca_result
    ))
}

plot_heterogeneity_pca <- function(full_input, pgs_vars, selected_pheno_vars, 
                                    mapper_pgs, mapper_pheno, plot_dir, train_idx, test_idx) {
    
    # Use the provided PCA model
    train_data <- full_input[train_idx, , drop = FALSE]
    test_data  <- full_input[test_idx, , drop = FALSE]
    pgs_train <- train_data[, pgs_vars]
    pgs_test <- test_data[, pgs_vars]
    pca_result <- prcomp(pgs_train, scale. = TRUE)
    pc_test <- predict(pca_result, newdata = pgs_test)
    # Variance explained plot on training set
    summ <- summary(pca_result)
    var_explained_prop <- summ$importance[2, ]
    df_var <- data.frame(
        PC = factor(paste0("PC", 1:12), levels = paste0("PC", 1:12)),
        Variance_Explained = var_explained_prop[1:12]
    )
    df_var$cumulative <- cumsum(df_var$Variance_Explained)
    p <- ggplot(df_var, aes(x = PC, y = Variance_Explained)) +
        geom_col(fill = "#d9642a", color = "#B22222") + 
        geom_text(aes(label = scales::percent(Variance_Explained, accuracy = 0.1)), 
                vjust = -0.5, size = 6, color = "#d9642a") +
        geom_line(aes(y = cumulative, group = 1), color = "#B22222", size = 1.2) +
        geom_point(aes(y = cumulative), color = "#B22222", size = 2) +
        geom_text(aes(y = cumulative, label = scales::percent(cumulative, accuracy = 0.1)), 
                vjust = -1.2, size = 6, color = "#B22222") +
        labs(
            x = "PGS principal components (PCs)",
            y = "Proportion of variance explained",
        ) +
        theme_bw(base_size = 18)
    ggsave(filename = file.path(plot_dir, "variance_explained_bar_plot_training.pdf"),
                plot = p, width = 8, height = 6, dpi = 100)
    # Make biplot colored by plaq_n
    # Rename loadings labels using mapper_pgs if provided
    rownames(pca_result$rotation) <- sapply(rownames(pca_result$rotation), mapper_pgs)
    loadings <- as.data.frame(pca_result$rotation[, 1:12])
    loadings_heatmap <- as.matrix(loadings)
    # Prepare matrix for ComplexHeatmap (PGS x PC)
    mat <- loadings_heatmap
    # Use hierarchical clustering with "ward.D2" method for better dendrograms
    # Take absolute value of loadings
    # Use ComplexHeatmap default clustering
    mat <- abs(mat)
    ht <- Heatmap(
        mat,
        name = "Loading",
        col = colorRamp2(
            c(0, max(mat)),
            c("white", "blue")
        ),
        row_names_side = "left",
        column_names_side = "bottom",
        column_names_rot = 90,
        column_title = "Principal component",
        row_title = "PGS",
        heatmap_legend_param = list(title = "Loading"),
        row_dend_side = "right", # Show row dendrogram on the right
        cluster_rows = readRDS(file.path(plot_dir, "pgs_effects_heatmap_stage_2_col_dend.rds")),
        cluster_columns = hclust(dist(t(mat)), method = "ward.D2")
    )
    pdf(file.path(plot_dir, "pgs_loadings_heatmap_predend.pdf"), width = 8, height = 4)
    draw(ht)
    dev.off()
    loadings$var <- rownames(loadings)
    # Add plaq_n to the pc_test data
    pc_test <- as.data.frame(pc_test)
    pc_test$plaq_n <- full_input[test_idx, "plaq_n"]
    pc_test <- pc_test[!is.na(pc_test$plaq_n), ]
    # Biplot for PC1 and PC2
    range_x12 <- max(abs(c(pc_test$PC1, loadings$PC1 * 7)))
    range_y12 <- max(abs(c(pc_test$PC2, loadings$PC2 * 7)))
    p_biplot_12 <- ggplot(pc_test, aes(x = PC1, y = PC2, color = plaq_n)) +
        geom_point(size = 3) +
        geom_segment(
            data = loadings,
            aes(x = 0, y = 0, xend = 7 * PC1, yend = 7 * PC2),
            arrow = arrow(length = unit(0.2, "cm")),
            color = "blue"
        ) +
        geom_text_repel(
            data = loadings,
            aes(x = 7 * PC1, y = 7 * PC2, label = var),
            size = 5, color = "black"
        ) +
        scale_color_gradient(low = "white", high = "red") +
        ggtitle("PCA Biplot of PGS Variables Colored by plaq_n (PC1 vs PC2)") +
        theme_bw(base_size = 16) +
        coord_cartesian(
            xlim = c(-range_x12, range_x12),
            ylim = c(-range_y12, range_y12),
            clip = "off"
        )

    ggsave(
        filename = file.path(plot_dir, "pca_biplot_colored_by_plaq_n_test_set_PC1_PC2.pdf"),
        plot = p_biplot_12, width = 8, height = 6, dpi = 100
    )

    # PC2 vs PC6
    range_x26 <- max(abs(c(pc_test$PC2, loadings$PC2 * 3)))
    range_y26 <- max(abs(c(pc_test$PC6, loadings$PC6 * 3)))

    p_biplot_26 <- ggplot(pc_test, aes(x = PC2, y = PC6, color = plaq_n)) +
        geom_point(size = 3) +
        geom_segment(
            data = loadings,
            aes(x = 0, y = 0, xend = 5 * PC2, yend = 5 * PC6),
            arrow = arrow(length = unit(0.2, "cm")),
            color = "blue"
        ) +
        geom_text_repel(
            data = loadings,
            aes(x = 5 * PC2, y = 5* PC6, label = var),
            size = 5, color = "black", max.overlaps=Inf
        ) +
        scale_color_gradient(low = "white", high = "red") +
        ggtitle("PCA Biplot of PGS Variables Colored by plaq_n (PC2 vs PC6)") +
        theme_bw(base_size = 16) +
        coord_cartesian(
            xlim = c(-range_x26, range_x26),
            ylim = c(-range_y26, range_y26),
            clip = "off"
        )
    ggsave(
        filename = file.path(plot_dir, "pca_biplot_colored_by_plaq_n_test_set_PC2_PC6.pdf"),
        plot = p_biplot_26, width = 8, height = 6, dpi = 100
    )
    range_x16 <- max(abs(c(pc_test$PC1, loadings$PC1 * 3)))
    range_y16 <- max(abs(c(pc_test$PC6, loadings$PC6 * 3)))
    p_biplot_16 <- ggplot(pc_test, aes(x = PC1, y = PC6, color = plaq_n)) +
        geom_point(size = 3) +
        geom_segment(
            data = loadings,
            aes(x = 0, y = 0, xend = 5 * PC1, yend = 5 * PC6),
            arrow = arrow(length = unit(0.2, "cm")),
            color = "blue"
        ) +
        geom_text_repel(
            data = loadings,
            aes(x = 5 * PC1, y = 5 * PC6, label = var),
            size = 5, color = "black", max.overlaps=Inf
        ) +
        scale_color_gradient(low = "white", high = "red") +
        ggtitle("PCA Biplot of PGS Variables Colored by plaq_n (PC1 vs PC6)") +
        theme_bw(base_size = 16) +
        coord_cartesian(
            xlim = c(-range_x16, range_x16),
            ylim = c(-range_y16, range_y16),
            clip = "off"
        )
    ggsave(
        filename = file.path(plot_dir, "pca_biplot_colored_by_plaq_n_test_set_PC1_PC6.pdf"),
        plot = p_biplot_16, width = 8, height = 6, dpi = 100
    )
}

plot_heterogeneity_pca_phenotypes <- function(full_input, selected_pheno_vars, 
                                   mapper_pgs, mapper_pheno, plot_dir) {
    
    # Impute missing values
    for (col in colnames(full_input)) {
        full_input[[col]][is.na(full_input[[col]])] <- mean(full_input[[col]], na.rm = TRUE)
    }
    phenotypes <- full_input[, selected_pheno_vars]
    
    pc_result <- prcomp(phenotypes, scale. = TRUE)

    # Variance explained plot on phenotypes
    summ <- summary(pc_result)
    var_explained_prop <- summ$importance[2, ]
    df_var <- data.frame(
        PC = factor(paste0("PC", 1:12), levels = paste0("PC", 1:12)),
        Variance_Explained = var_explained_prop[1:12]
    )
    p <- ggplot(df_var, aes(x = PC, y = Variance_Explained)) +
        geom_col(color="blue") + 
        geom_text(aes(label = scales::percent(Variance_Explained, accuracy = 0.1)), 
                          vjust = -0.5, size = 5) +
        labs(
            x = "Principal Components",
            y = "Proportion of Variance Explained"
        ) +
        theme_bw(base_size = 18)
    ggsave(filename = file.path(plot_dir, "variance_explained_bar_plot_phenotypes.pdf"),
             plot = p, width = 8, height = 6, dpi = 100)

    # PCA on phenotypes was just computed as pc_result
    # Color by FH1263
    color_var <- "FH1263_SCORE1_AVG"

    # Extract scores and attach the coloring variable
    scores <- as.data.frame(pc_result$x[, 1:2])
    colnames(scores) <- c("PC1", "PC2")
    scores$FH1263 <- full_input[, color_var]

    # Extract loadings for the first two PCs
    loadings <- as.data.frame(pc_result$rotation[, 1:2])
    colnames(loadings) <- c("PC1", "PC2")
    selected_phenotype_variables <- c(
        "amyloid", "nft", "plaq_n", "plaq_d", "tangles", "gpath", "cogn_global_lv"
    )
    # zero out PC1 and PC2 for any loading not in selected_phenotype_variables
    loadings[!rownames(loadings) %in% selected_phenotype_variables,
             c("PC1", "PC2")] <- 0
    loadings$var <- ifelse(
        rownames(loadings) %in% selected_phenotype_variables,
        rownames(loadings),
        ""
    )
    loadings$var <- sapply(loadings$var, mapper_pheno)

    # Determine limits for arrows and points
    scale_factor <- 35
    x_lim <- max(abs(c(scores$PC1, loadings$PC1 * scale_factor)))
    y_lim <- max(abs(c(scores$PC2, loadings$PC2 * scale_factor)))
    # Z-transform FH1263
    scores$FH1263 <- scale(scores$FH1263)
    # Build the biplot
    p_biplot_pheno <- ggplot(scores, aes(x = PC1, y = PC2, color = FH1263)) +
        geom_point(size = 1) +
        geom_segment(
            data = loadings,
            aes(x = 0, y = 0, xend = PC1 * scale_factor, yend = PC2 * scale_factor),
            arrow = arrow(length = unit(0.2, "cm")),
            color = "red"
        ) +
        geom_text_repel(
            data = loadings,
            aes(x = PC1 * scale_factor, y = PC2 * scale_factor, label = var),
            size = 5, color = "black", max.overlaps=Inf
        ) +
        scale_color_gradient(low = "white", high = "blue") +
        coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(-y_lim, y_lim), clip = "off") +
        theme_bw(base_size = 16)

    # Save to file
    ggsave(
        filename = file.path(plot_dir, glue("pca_biplot_phenotypes_colored_by_{color_var}_final.pdf")),
        plot = p_biplot_pheno,
        width = 8,
        height = 6,
        dpi = 100
    )
  
}

cluster_individuals_by_pgs <- function(full_input, pgs_vars, mapper_pgs, plot_dir, train_idx, test_idx) {
    train_data <- full_input[train_idx, , drop = FALSE]
    test_data  <- full_input[test_idx, , drop = FALSE]
    pgs_train <- train_data[, pgs_vars]
    pca_result <- prcomp(pgs_train, scale. = TRUE)
    summ <- summary(pca_result)
    var_explained_prop <- summ$importance[2, ]
    df_var <- data.frame(
        PC = factor(paste0("PC", 1:12), levels = paste0("PC", 1:12)),
        Variance_Explained = var_explained_prop[1:12]
    )
    pgs_data <- full_input[, pgs_vars, drop = FALSE]
    pgs_scaled <- scale(pgs_data)
    pc_scores <- predict(pca_result, newdata = pgs_data)
    pc_scaled <- sweep(pc_scores, 2, var_explained_prop[1:12], FUN = "*")
    wss <- sapply(1:10, function(k) {
        kmeans(pc_scaled, centers = k, nstart = 10)$tot.withinss
    })
    elbow_df <- data.frame(
        k = 1:10,
        wss = wss
    )
    p_elbow <- ggplot(elbow_df, aes(x = k, y = wss)) +
        geom_line(color = "blue") +
        geom_point(color = "red", size = 2) +
        scale_x_continuous(breaks = 1:10) +
        labs(
            x = "Number of Clusters (k)",
            y = "Total Within-Cluster Sum of Squares (WSS)",
            title = "Elbow Method for Determining Optimal k"
        ) +
        theme_bw(base_size = 18)
    ggsave(
        filename = file.path(plot_dir, "elbow_method_kmeans_pgs.pdf"),
        plot = p_elbow,
        width = 8,
        height = 6,
        dpi = 100
    )
set.seed(123)
k_opt <- 6
kmeans_result <- kmeans(pc_scaled, centers = k_opt, nstart = 25)
clusters <- kmeans_result$cluster
# Project all data onto PCA space
pc_all <- pc_scores
pc_df <- as.data.frame(pc_all[, 1:3])
pc_df$cluster <- factor(clusters)

# Get loadings for biplot
loadings <- as.data.frame(pca_result$rotation[, 1:3])
rownames(loadings) <- sapply(rownames(loadings), mapper_pgs)
loadings$var <- rownames(loadings)

# Create biplot with clusters (PC1 vs PC2)
range_x12 <- max(abs(c(pc_df$PC1, loadings$PC1 * 7)))
range_y12 <- max(abs(c(pc_df$PC2, loadings$PC2 * 7)))
p_biplot_clusters_12 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_segment(
        data = loadings,
        aes(x = 0, y = 0, xend = 7 * PC1, yend = 7 * PC2),
        arrow = arrow(length = unit(0.2, "cm")),
        color = "black",
        inherit.aes = FALSE
    ) +
    geom_text_repel(
        data = loadings,
        aes(x = 7 * PC1, y = 7 * PC2, label = var),
        size = 6, color = "black",
        inherit.aes = FALSE
    ) +
    scale_color_discrete(name = "Cluster") +
    labs(
        title = "PCA Biplot with K-means Clusters (k=6, PC1 vs PC2)",
        x = "PC1 (26.9%)",
        y = "PC2 (15.4%)"
    ) +
    theme_bw(base_size = 32) +
    theme(
        text = element_text(color = "black", size = 32),
        axis.text = element_text(color = "black", size = 32),
        axis.title = element_text(color = "black", size = 32),
        plot.title = element_text(color = "black", size = 32),
        legend.text = element_text(color = "black", size = 32),
        legend.title = element_text(color = "black", size = 32)
    ) +
    coord_cartesian(
        xlim = c(-range_x12, range_x12),
        ylim = c(-range_y12, range_y12),
        clip = "off"
    )

ggsave(
    filename = file.path(plot_dir, "pca_biplot_kmeans_clusters_PC1_PC2.pdf"),
    plot = p_biplot_clusters_12,
    width = 10,
    height = 8,
    dpi = 300
)

# Save cluster assignments
cluster_results <- data.frame(
    individual_id = full_input$X.projid,
    cluster = clusters
)
write.table(
    cluster_results,
    file = file.path(plot_dir, "kmeans_cluster_assignments.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
# Read in cluster assignments from previous analysis
cluster_results <- read.table(
    file.path(plot_dir, "kmeans_cluster_assignments.tsv"),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
)

# Extract cluster assignments
clusters <- cluster_results$cluster

# Rank clusters by average global AD pathology (gpath)
gpath_means <- tapply(full_input$gpath, clusters, mean, na.rm = TRUE)
gpath_ranking <- order(gpath_means, decreasing = TRUE)
gpath_ranked_means <- gpath_means[gpath_ranking]

# Print ranking
cat("Clusters ranked by average global AD pathology (gpath):\n")
for (i in 1:length(gpath_ranked_means)) {
    cluster_id <- names(gpath_ranked_means)[i]
    mean_val <- round(gpath_ranked_means[i], 3)
    cat("Rank", i, ": Cluster", cluster_id, "- Mean gpath:", mean_val, "\n")
    cat("Number of individuals in Cluster", cluster_id, ":", sum(clusters == cluster_id), "\n")
}
k_opt <- length(unique(clusters))
# Create cross-tabulation of clusters and ApoE status
cluster_apoe_table <- table(clusters, full_input$apoe_genotype)
cluster_apoe_df <- as.data.frame(cluster_apoe_table)
colnames(cluster_apoe_df) <- c("Cluster", "ApoE_Genotype", "Count")
# Reorder clusters by average gpath means
cluster_order <- names(sort(gpath_means, decreasing = TRUE))
cluster_apoe_df$Cluster <- factor(cluster_apoe_df$Cluster, levels = cluster_order)

# Create stacked bar plot
p_cluster_apoe <- ggplot(cluster_apoe_df, aes(x = Cluster, y = Count, fill = factor(ApoE_Genotype))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
        x = "Cluster (ordered by average gpath)",
        y = "Count",
        fill = "ApoE Genotype",
        title = "Cluster Assignments by ApoE Genotype"
    ) +
    theme_bw(base_size = 32) +
    scale_fill_brewer(type = "qual", palette = "Set2", 
                      labels = c("e2e2", "e2e3", "e2e4", "e3e3", "e3e4", "e4e4"))

ggsave(
    filename = file.path(plot_dir, "cluster_apoe_concordance_stacked.pdf"),
    plot = p_cluster_apoe,
    width = 8,
    height = 6,
    dpi = 300
)

# Calculate cluster-level averages for key phenotypes
phenotype_vars <- c("cogn_global_lv", "gpath", "plaq_n", "amyloid", "plaq_d", "tangles")
phenotype_labels <- c("Global cognitive function (19 tests)", "Global AD pathology burden", "Neuritic plaque burden (5 regions)", "Amyloid level (% cortex area, 8 brain regions)", "Diffuse plaque burden (5 regions)", "Tangle density (IHC, 8 brain regions)")

cluster_phenotype_means <- data.frame(
    Cluster = 1:k_opt,
    stringsAsFactors = FALSE
)

# Z-score the phenotype variables before calculating means
# Create a copy of full_input with z-scored phenotype variables
full_input_zscore <- full_input
for (var_name in phenotype_vars) {
    full_input_zscore[[var_name]] <- scale(full_input[[var_name]])[,1]
}

# Calculate means using the z-scored data
for (i in 1:length(phenotype_vars)) {
    var_name <- phenotype_vars[i]
    means_by_cluster <- tapply(full_input_zscore[[var_name]], clusters, mean, na.rm = TRUE)
    cluster_phenotype_means[[phenotype_labels[i]]] <- means_by_cluster
}

# Convert to matrix for heatmap
rownames(cluster_phenotype_means) <- paste0("Cluster ", cluster_phenotype_means$Cluster)
heatmap_matrix <- as.matrix(cluster_phenotype_means[, -1])

# Reorder clusters by gpath means (already calculated above)
cluster_order_by_gpath <- order(gpath_means, decreasing = TRUE)
heatmap_matrix_reordered <- heatmap_matrix[cluster_order_by_gpath, ]

# Create heatmap with clusters reordered by gpath
ht_cluster_phenotypes <- Heatmap(
    heatmap_matrix_reordered,
    name = "Cluster average of standardized ROSMAP AD phenotype",
    col = colorRamp2(c(-.4, 0, .4), c("violet", "white", "red")),
    row_names_side = "left",
    column_names_side = "bottom",
    column_names_rot = 90,
    column_title = "Selected ROSMAP AD phenotype",
    row_title = "AD genetic liability profile subtype",
    heatmap_legend_param = list(title = "Cluster average of standardized ROSMAP AD phenotype"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_max_width = unit(4, "cm")
)

pdf(file.path(plot_dir, "cluster_phenotype_averages_heatmap.pdf"), width = 8, height = 6)
draw(ht_cluster_phenotypes)
dev.off()


# Create violin plots comparing cluster 5 vs cluster 1 for neuritic and diffuse plaque burden
clusters_to_compare <- c(2,3)
plaque_vars <- c("plaq_n", "plaq_d")
plaque_labels <- c("Neuritic plaque burden", "Diffuse plaque burden")

# Create data for violin plots
violin_data <- data.frame(
    individual_id = full_input_zscore$X.projid,
    cluster = clusters,
    plaq_n = full_input_zscore$plaq_n,
    plaq_d = full_input_zscore$plaq_d
)

# Filter for clusters 1 and 5 only
violin_data_filtered <- violin_data[violin_data$cluster %in% clusters_to_compare, ]
violin_data_filtered$cluster <- factor(violin_data_filtered$cluster)
# Perform statistical tests
test_plaq_n <- wilcox.test(plaq_n ~ cluster, data = violin_data_filtered)
test_plaq_d <- wilcox.test(plaq_d ~ cluster, data = violin_data_filtered)

# Print p-values
cat("Neuritic plaque burden (plaq_n) - Wilcoxon test p-value:", test_plaq_n$p.value, "\n")
cat("Diffuse plaque burden (plaq_d) - Wilcoxon test p-value:", test_plaq_d$p.value, "\n")
# Extract data for clusters 2 and 3
cluster2_data <- violin_data_filtered[violin_data_filtered$cluster == 2, ]
cluster3_data <- violin_data_filtered[violin_data_filtered$cluster == 3, ]

# Compute and print all 4 medians
cat("Cluster 2 neuritic plaque burden median:", median(cluster2_data$plaq_n, na.rm = TRUE), "\n")
cat("Cluster 2 diffuse plaque burden median:", median(cluster2_data$plaq_d, na.rm = TRUE), "\n")
cat("Cluster 3 neuritic plaque burden median:", median(cluster3_data$plaq_n, na.rm = TRUE), "\n")
cat("Cluster 3 diffuse plaque burden median:", median(cluster3_data$plaq_d, na.rm = TRUE), "\n")

# Alternative version: Group by plaque type
# Create data for alternative violin plot (plaque type grouped together)
cluster2_neuritic <- data.frame(
    value = cluster2_data$plaq_n,
    cluster = "Cluster 2",
    type = "Neuritic"
)
cluster3_neuritic <- data.frame(
    value = cluster3_data$plaq_n,
    cluster = "Cluster 3", 
    type = "Neuritic"
)
cluster2_diffuse <- data.frame(
    value = cluster2_data$plaq_d,
    cluster = "Cluster 2",
    type = "Diffuse"
)
cluster3_diffuse <- data.frame(
    value = cluster3_data$plaq_d,
    cluster = "Cluster 3",
    type = "Diffuse"
)

grouped_data <- rbind(cluster2_neuritic, cluster3_neuritic, cluster2_diffuse, cluster3_diffuse)
grouped_data$cluster <- factor(grouped_data$cluster)
grouped_data$type <- factor(grouped_data$type)

p_violin_grouped <- ggplot(grouped_data, aes(x = type, y = value, fill = cluster)) +
    geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(width = 0.8)) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), alpha = 0.8) +
    stat_summary(fun = "mean", geom = "point", size = 3, color = "black", 
                 position = position_dodge(width = 0.8)) +
    labs(
        x = "Plaque Type",
        y = "Plaque Burden",
        title = "Neuritic vs Diffuse Plaque Burden by Cluster",
        fill = "Cluster"
    ) +
    scale_fill_discrete(type = scales::hue_pal()(6)[c(2, 3)]) +
    ylim(-2, 4) +
    theme_bw(base_size = 16)

ggsave(
    filename = file.path(plot_dir, "violin_grouped_plaque_by_type.pdf"),
    plot = p_violin_grouped,
    width = 8,
    height = 6,
    dpi = 300
)

# Stratify by APOE genotype
apoe_genotypes <- unique(full_input$apoe_genotype)

# For each APOE genotype, create the violin plot analysis
for (apoe_geno in apoe_genotypes) {
    # Filter data for current APOE genotype
    apoe_subset_idx <- which(full_input$apoe_genotype == apoe_geno)
    violin_data_apoe <- violin_data[violin_data$individual_id %in% full_input$X.projid[apoe_subset_idx], ]
    
    # Filter for clusters 2 and 3 only
    violin_data_apoe_filtered <- violin_data_apoe[violin_data_apoe$cluster %in% clusters_to_compare, ]
    
    # Skip if not enough data
    if (nrow(violin_data_apoe_filtered) < 5) {
        cat("Skipping APOE genotype", apoe_geno, "- insufficient data\n")
        next
    }
    
    violin_data_apoe_filtered$cluster <- factor(violin_data_apoe_filtered$cluster)
    
    # Perform statistical tests for this APOE genotype
    test_plaq_n_apoe <- wilcox.test(plaq_n ~ cluster, data = violin_data_apoe_filtered)
    test_plaq_d_apoe <- wilcox.test(plaq_d ~ cluster, data = violin_data_apoe_filtered)
    
    cat("\nAPOE genotype:", apoe_geno, "\n")
    cat("Neuritic plaque burden (plaq_n) - Wilcoxon test p-value:", test_plaq_n_apoe$p.value, "\n")
    cat("Diffuse plaque burden (plaq_d) - Wilcoxon test p-value:", test_plaq_d_apoe$p.value, "\n")
    
    # Extract data for clusters 2 and 3 within this APOE genotype
    cluster2_data_apoe <- violin_data_apoe_filtered[violin_data_apoe_filtered$cluster == 2, ]
    cluster3_data_apoe <- violin_data_apoe_filtered[violin_data_apoe_filtered$cluster == 3, ]
    
    # Compute and print medians
    cat("Cluster 2 neuritic plaque burden median:", median(cluster2_data_apoe$plaq_n, na.rm = TRUE), "\n")
    cat("Cluster 2 diffuse plaque burden median:", median(cluster2_data_apoe$plaq_d, na.rm = TRUE), "\n")
    cat("Cluster 3 neuritic plaque burden median:", median(cluster3_data_apoe$plaq_n, na.rm = TRUE), "\n")
    cat("Cluster 3 diffuse plaque burden median:", median(cluster3_data_apoe$plaq_d, na.rm = TRUE), "\n")
    cat("Number of individuals in Cluster 2:", nrow(cluster2_data_apoe), "\n")
    cat("Number of individuals in Cluster 3:", nrow(cluster3_data_apoe), "\n")
    
    # Create grouped data for violin plot
    cluster2_neuritic_apoe <- data.frame(
        value = cluster2_data_apoe$plaq_n,
        cluster = "Cluster 2",
        type = "Neuritic"
    )
    cluster3_neuritic_apoe <- data.frame(
        value = cluster3_data_apoe$plaq_n,
        cluster = "Cluster 3", 
        type = "Neuritic"
    )
    cluster2_diffuse_apoe <- data.frame(
        value = cluster2_data_apoe$plaq_d,
        cluster = "Cluster 2",
        type = "Diffuse"
    )
    cluster3_diffuse_apoe <- data.frame(
        value = cluster3_data_apoe$plaq_d,
        cluster = "Cluster 3",
        type = "Diffuse"
    )
    
    grouped_data_apoe <- rbind(cluster2_neuritic_apoe, cluster3_neuritic_apoe, 
                               cluster2_diffuse_apoe, cluster3_diffuse_apoe)
    grouped_data_apoe$cluster <- factor(grouped_data_apoe$cluster)
    grouped_data_apoe$type <- factor(grouped_data_apoe$type)
    
    p_violin_grouped_apoe <- ggplot(grouped_data_apoe, aes(x = type, y = value, fill = cluster)) +
        geom_violin(trim = FALSE, alpha = 0.7, position = position_dodge(width = 0.8)) +
        geom_boxplot(width = 0.2, position = position_dodge(width = 0.8), alpha = 0.8) +
        stat_summary(fun = "mean", geom = "point", size = 3, color = "black", 
                     position = position_dodge(width = 0.8)) +
        labs(
            x = "Plaque Type",
            y = "Plaque Burden",
            title = glue("Neuritic vs Diffuse Plaque Burden by Cluster (APOE {apoe_geno})"),
            fill = "Cluster"
        ) +
        scale_fill_discrete(type = scales::hue_pal()(6)[c(2, 3)]) +
        ylim(-2, 4) +
        theme_bw(base_size = 16)
    
    ggsave(
        filename = file.path(plot_dir, glue("violin_grouped_plaque_by_type_apoe_{apoe_geno}.pdf")),
        plot = p_violin_grouped_apoe,
        width = 8,
        height = 6,
        dpi = 300
    )
}
}