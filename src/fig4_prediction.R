# src/fig4_prediction.R -- Figure 4
# Nested 5-fold CV with a partially penalised ridge: covariates and APOE unpenalised, the PGS
# block shrunk with a penalty chosen by inner 5-fold CV. Covariate imputation means and the PGS
# PCA are refit inside every outer fold; outcomes are without imputation.
#   run_cv5()                  out-of-fold predictions, pooled and per-fold r, paired bootstrap vs APOE
#   make_figure4_prediction()  Fig4_delta_r_vs_apoe.pdf and Table S10
suppressPackageStartupMessages(library(ggplot2))

CV_K <- 5
CV_B <- 4000
CV_LAMBDA_GRID <- 10^seq(-2, 6, length.out = 81)
CV_FH_PGS <- 'FH1263_SCORE1_AVG'
CV_MODEL_COLOURS <- c('Covariates only' = '#B09C85', 'APOE' = '#7E6148', 'APOE + AD PGS' = '#3C5488',
                      'APOE + All-PGS' = '#00A087', 'APOE + PGS PCs' = '#91D1C2')
CV_MODEL_DISPLAY <- list('Covariates only', expression(italic(APOE)), expression(italic(APOE) * ' + AD PGS'),
                         expression(italic(APOE) * ' + All-PGS'), expression(italic(APOE) * ' + PGS PCs'))
CV_PHENOTYPE_SHORT <- c(cogn_global_lv = 'Global\ncognition', gpath = 'Global\npathology', plaq_n = 'Neuritic\nplaques',
                        amyloid = 'Amyloid', plaq_d = 'Diffuse\nplaques', tangles = 'Tangle\ndensity', nft = 'NFT\nburden',
                        niareagansc = 'NIA-Reagan\nscore', cogdx = 'Cognitive\ndiagnosis')
# Exact partially penalised ridge: lambda added only to the training-standardised PGS coordinates.
partial_ridge_fit <- function(data, model, train_rows, y) {
    X_unpenalized <- cbind(intercept = 1, as.matrix(data[, model$unpenalized, drop = FALSE]))
    if (length(model$penalized) == 0) {
        beta <- solve(crossprod(X_unpenalized[train_rows, , drop = FALSE]),
                      crossprod(X_unpenalized[train_rows, , drop = FALSE], y[train_rows]))
        return(function(lambda, rows_new) drop(X_unpenalized[rows_new, , drop = FALSE] %*% beta))
    }
    X_pen_raw <- as.matrix(data[, model$penalized, drop = FALSE])
    centres <- colMeans(X_pen_raw[train_rows, , drop = FALSE])
    scales <- apply(X_pen_raw[train_rows, , drop = FALSE], 2, sd)
    X <- cbind(X_unpenalized, scale(X_pen_raw, centres, scales))
    XtX <- crossprod(X[train_rows, , drop = FALSE])
    Xty <- crossprod(X[train_rows, , drop = FALSE], y[train_rows])
    penalized_index <- ncol(X_unpenalized) + seq_along(model$penalized)
    function(lambda, rows_new) {
        A <- XtX
        diag(A)[penalized_index] <- diag(A)[penalized_index] + lambda
        drop(X[rows_new, , drop = FALSE] %*% solve(A, Xty))
    }
}

cv_theme <- function() {
    theme_house() + theme(legend.position = 'top', panel.grid.major.x = element_blank(),
                          plot.caption = element_text(size = HOUSE_FONT_PT, colour = HOUSE_TEXT_COLOUR))
}

# Nested 5-fold CV. One row per model x phenotype: pooled and per-fold r, paired bootstrap vs APOE.
run_cv5 <- function(seed = SPLIT_SEED, B = CV_B) {
    t0 <- Sys.time()
    K <- CV_K
    cohort <- load_cohort()
    pgs_columns <- grep('_SCORE1_AVG$', names(cohort), value = TRUE)
    stopifnot(length(pgs_columns) == 12, CV_FH_PGS %in% pgs_columns)
    pc_columns <- paste0('PGSPC', 1:N_PGS_PCS)
    feature_columns <- c(pgs_columns, COVARIATES, 'apoe_genotype')
    n_individuals <- nrow(cohort)

    MODELS <- list(
        'Covariates only' = list(unpenalized = COVARIATES, penalized = character(0)),
        'APOE' = list(unpenalized = c('apoe_genotype', COVARIATES), penalized = character(0)),
        'APOE + AD PGS' = list(unpenalized = c('apoe_genotype', COVARIATES), penalized = CV_FH_PGS),
        'APOE + All-PGS' = list(unpenalized = c('apoe_genotype', COVARIATES), penalized = pgs_columns),
        'APOE + PGS PCs' = list(unpenalized = c('apoe_genotype', COVARIATES), penalized = pc_columns))
    BASELINE <- 'APOE'
    COMPARED <- names(MODELS)[3:5]

    # ---- Out-of-fold predictions: everything refit without fold k ----------------------
    set.seed(seed)
    outer_fold <- sample(rep(seq_len(K), length.out = n_individuals))
    oof <- lapply(setNames(nm = PREDICTION_PHENOTYPES), function(phenotype) {
        predictions <- matrix(NA_real_, n_individuals, length(MODELS), dimnames = list(NULL, names(MODELS)))
        for (k in seq_len(K)) {
            data <- cohort
            in_train <- outer_fold != k
            means <- colMeans(data[in_train, feature_columns], na.rm = TRUE)
            for (column in feature_columns) data[[column]][is.na(data[[column]])] <- means[column]
            pca <- prcomp(data[in_train, pgs_columns], scale. = TRUE)
            data[pc_columns] <- predict(pca, data[, pgs_columns])[, 1:N_PGS_PCS]

            observed <- !is.na(data[[phenotype]])
            train_rows <- which(in_train & observed)
            test_rows <- which(!in_train & observed)
            y <- data[[phenotype]]

            set.seed(seed * 100 + k)
            inner_fold <- sample(rep(seq_len(K), length.out = length(train_rows)))
            for (model_name in names(MODELS)) {
                model <- MODELS[[model_name]]
                if (length(model$penalized) == 0) {
                    lambda <- NA_real_
                } else {
                    inner_r <- rowMeans(sapply(seq_len(K), function(j) {
                        fit <- partial_ridge_fit(data, model, train_rows[inner_fold != j], y)
                        vapply(CV_LAMBDA_GRID, function(l)
                            cor(fit(l, train_rows[inner_fold == j]), y[train_rows[inner_fold == j]]), numeric(1))
                    }))
                    lambda <- CV_LAMBDA_GRID[which.max(inner_r)]
                }
                fit <- partial_ridge_fit(data, model, train_rows, y)
                predictions[test_rows, model_name] <- fit(lambda, test_rows)
            }
        }
        predictions
    })
    say('CV: out-of-fold predictions done ({round(as.numeric(Sys.time() - t0, units = "secs"))} s)')

    # ---- Pooled r, per-fold r, paired bootstrap contrasts vs APOE ------------------------
    do.call(rbind, lapply(PREDICTION_PHENOTYPES, function(phenotype) {
        y_all <- cohort[[phenotype]]
        observed <- which(!is.na(y_all))
        predictions <- oof[[phenotype]][observed, ]
        y <- y_all[observed]
        fold <- outer_fold[observed]
        pooled_r <- cor(predictions, y)[, 1]
        fold_r <- sapply(seq_len(K), function(k) cor(predictions[fold == k, ], y[fold == k])[, 1])

        set.seed(seed * 7 + match(phenotype, PREDICTION_PHENOTYPES))
        boot_r <- replicate(B, {
            i <- sample(length(y), replace = TRUE)
            cor(predictions[i, ], y[i])[, 1]
        })
        delta_boot <- boot_r[COMPARED, , drop = FALSE] - matrix(boot_r[BASELINE, ], length(COMPARED), B, byrow = TRUE)
        p_one <- rowMeans(delta_boot > 0)

        per_model <- data.frame(phenotype = phenotype, model = names(MODELS), n_observed = length(observed),
                                pooled_oof_r = pooled_r, fold_r, row.names = NULL)
        names(per_model)[5:(4 + K)] <- paste0('fold', seq_len(K))
        per_model$delta_vs_APOE <- NA; per_model$se_delta <- NA
        per_model$ci_low <- NA; per_model$ci_high <- NA; per_model$p_two_sided <- NA
        idx <- match(COMPARED, per_model$model)
        per_model$delta_vs_APOE[idx] <- pooled_r[COMPARED] - pooled_r[[BASELINE]]
        per_model$se_delta[idx] <- apply(delta_boot, 1, sd)
        per_model$ci_low[idx] <- apply(delta_boot, 1, quantile, 0.025)
        per_model$ci_high[idx] <- apply(delta_boot, 1, quantile, 0.975)
        per_model$p_two_sided[idx] <- pmax(2 * pmin(p_one, 1 - p_one), 1 / B)
        per_model
    }))
}

make_figure4_prediction <- function() {
    t0 <- Sys.time()
    MODEL_NAMES <- names(CV_MODEL_COLOURS)
    BASELINE <- 'APOE'
    COMPARED <- MODEL_NAMES[3:5]
    results <- run_cv5(SPLIT_SEED, CV_B)

    results$phenotype_label <- factor(CV_PHENOTYPE_SHORT[results$phenotype], levels = CV_PHENOTYPE_SHORT)
    results$model <- factor(results$model, levels = MODEL_NAMES)

    # ---- Contrasts: bar = pooled delta r vs APOE, whisker = paired 95% bootstrap CI ---------------
    delta_rows <- subset(results, model %in% COMPARED)
    delta_rows$stars <- ifelse(delta_rows$p_two_sided < 0.001, '***', ifelse(delta_rows$p_two_sided < 0.01, '**',
                        ifelse(delta_rows$p_two_sided < 0.05, '*', '')))
    delta_rows$star_y <- ifelse(delta_rows$delta_vs_APOE >= 0, delta_rows$ci_high + 0.002, delta_rows$ci_low - 0.002)
    delta_rows$star_hjust <- ifelse(delta_rows$delta_vs_APOE >= 0, 0, 1)
    deltas <- ggplot(delta_rows, aes(x = phenotype_label, y = delta_vs_APOE, fill = model)) +
        geom_hline(yintercept = 0, linewidth = HOUSE_LINE_MM, colour = 'grey40') +
        geom_col(position = position_dodge(width = 0.8), width = 0.7) +
        geom_errorbar(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.8),
                      width = 0.25, linewidth = HOUSE_LINE_MM, colour = 'grey25') +
        geom_text(aes(y = star_y, label = stars, hjust = star_hjust), position = position_dodge(width = 0.8),
                  angle = 90, size = HOUSE_FONT_PT / .pt, family = HOUSE_FONT_FAMILY, colour = HOUSE_TEXT_COLOUR) +
        scale_fill_manual(values = CV_MODEL_COLOURS[COMPARED], name = NULL, labels = CV_MODEL_DISPLAY[3:5]) +
        guides(fill = guide_legend(nrow = 1)) +
        labs(x = 'Observed ROSMAP Phenotypes', y = expression(Delta * italic(r) * ' vs ' * italic(APOE)),
             caption = '* p<0.05, ** p<0.01, *** p<0.001 (two-sided paired bootstrap)') +
        cv_theme()
    deltas_path <- fig_path('Figure4', 'Fig4_delta_r_vs_apoe.pdf')
    open_panel_pdf(deltas_path, 'prediction'); print(deltas); invisible(dev.off())

    # ---- Table S10: the long results table with display names -----------------------------------
    table_s10 <- results[, c('phenotype', 'model', 'n_observed', 'pooled_oof_r', paste0('fold', 1:CV_K),
                             'delta_vs_APOE', 'se_delta', 'p_two_sided')]
    table_s10$phenotype <- unname(mapper_pheno(as.character(table_s10$phenotype)))
    names(table_s10) <- c('Observed ROSMAP Phenotype', 'Model Predictors', 'n_observed',
                          'Pooled out-of-fold test split correlation (r)',
                          paste0('Fold', 1:CV_K, ' test split correlation (r)'),
                          'Delta r vs APOE', 'Delta r standard error (bootstrap)', 'P-value (two-sided, bootstrap)')
    s10_path <- table_path('TableS10_predictive_models.tsv')
    write.table(table_s10, s10_path, sep = '\t', row.names = FALSE, quote = FALSE)

    say('Figure 4 / Table S10 done ({round(as.numeric(Sys.time() - t0, units = "secs"))} s)')
    rbind(
        manifest_row('Figure4', 'delta r vs APOE', deltas_path, note = glue('nested 5-fold CV (seed {SPLIT_SEED}), paired bootstrap 95% CI (B = {CV_B})')),
        manifest_row('TableS10', 'table', s10_path, note = 'per model x phenotype: n, pooled and per-fold out-of-fold r, delta r vs APOE, bootstrap SE and p'))
}
