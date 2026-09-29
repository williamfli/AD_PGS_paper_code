# src/data.R -- shared loaders: label mappers, the association result matrices, the cohort
# table, and the seeded split + PGS PCA that Figure 5, S5 and S14 share. Results are memoised.

.cache <- new.env(parent = emptyenv())
memo <- function(key, expr) {
    if (!exists(key, envir = .cache, inherits = FALSE)) assign(key, expr, envir = .cache)
    get(key, envir = .cache, inherits = FALSE)
}

# ---- Label mappers (ids -> display names) ------------------------------------------------
make_label_mapper <- function(annotation_df) {
    mapping <- setNames(annotation_df$name, annotation_df$variable)
    function(variable) {
        if (length(variable) != 1) return(vapply(variable, function(v) if (v %in% names(mapping)) mapping[[v]] else v, character(1)))
        if (variable %in% names(mapping)) mapping[[variable]] else variable
    }
}

# phenotype id (plaq_n) -> "Neuritic plaque burden (5 regions)"
mapper_pheno <- make_label_mapper(read.delim(input_path('phenotype_annotation'), header = TRUE, sep = '\t', fill = TRUE))

mapper_pgs <- local({
    ann <- read.delim(input_path('pgs_annotation'), header = TRUE, sep = '\t', fill = TRUE)
    names(ann) <- as.character(unlist(ann[1, ]))
    ann <- ann[-1, ]
    names(ann)[names(ann) == 'Trait Name'] <- 'name'
    names(ann)[names(ann) == 'GBE ID'] <- 'variable'
    ann$variable <- paste0(ann$variable, '_SCORE1_AVG')
    make_label_mapper(ann)
})
mapper_pgs_id <- function(id) mapper_pgs(paste0(id, '_SCORE1_AVG'))

# ---- Association result matrices (rows = PGS, columns = phenotypes) -------------------------
# Written by make_associations() (src/fig2_associations.R) into OUTPUT_DIR/association/<run>/.

ASSOCIATION_RUNS <- c('stage1', 'stage2', 'stage2_no_apoe', 'apoe_only')
result_dir <- function(run) {
    stopifnot(run %in% ASSOCIATION_RUNS)
    dir <- glue('{OUTPUT_DIR}association/{run}/')
    if (!all(file.exists(glue('{dir}{c("effects", "errors", "pvals", "fdr")}.tsv'))))
        stop(glue('association run "{run}" not found under {dir}; run `Rscript scripts/run_all.R Associations` first'))
    dir
}

read_result_dir <- function(run, mapped = TRUE) {
    dir <- result_dir(run)
    read_one <- function(name) {
        m <- read.delim(glue('{dir}{name}.tsv'), header = TRUE, sep = '\t', fill = TRUE)
        if (mapped) {
            rownames(m) <- sapply(rownames(m), mapper_pgs)
            colnames(m) <- sapply(colnames(m), mapper_pheno)
        }
        m
    }
    list(effects = read_one('effects'), errors = read_one('errors'), fdr = read_one('fdr'), pvals = read_one('pvals'))
}

load_stage1         <- function(mapped = TRUE) memo(glue('stage1_{mapped}'),   read_result_dir('stage1', mapped))
load_stage2         <- function(mapped = TRUE) memo(glue('stage2_{mapped}'),   read_result_dir('stage2', mapped))
load_stage2_no_apoe <- function(mapped = TRUE) memo(glue('stage2na_{mapped}'), read_result_dir('stage2_no_apoe', mapped))

# Stage 1 restricted to the PGS with any FDR < 0.5 (the 12 prioritized PGS; Fig 2a).
load_stage1_prioritized <- function(mapped = TRUE) {
    s1 <- load_stage1(mapped)
    keep <- rownames(s1$fdr)[apply(s1$fdr, 1, function(row) any(row < FDR_STAGE1, na.rm = TRUE))]
    lapply(s1, function(m) m[keep, , drop = FALSE])
}

load_apoe_only <- function(mapped = TRUE) memo(glue('apoeonly_{mapped}'), {
    read_one <- function(name) {
        m <- read.delim(glue('{result_dir("apoe_only")}{name}.tsv'), header = TRUE, sep = '\t', fill = TRUE)
        if (mapped) colnames(m) <- sapply(colnames(m), mapper_pheno)
        m
    }
    list(effects = read_one('effects'), errors = read_one('errors'), fdr = read_one('fdr'))
})

pgs_variable_ids     <- function() memo('pgs_ids',   rownames(load_stage2(mapped = FALSE)$effects))
stage2_phenotype_ids <- function() memo('pheno_ids', colnames(load_stage2(mapped = FALSE)$effects))

# ---- The cohort table (one row per participant, in the paper's row order) -------------------
load_cohort <- function() memo('cohort', {
    d <- read.delim(input_path('cohort_table'), header = TRUE, sep = '\t', fill = TRUE, colClasses = c('X.projid' = 'character'))
    d$X.projid <- formatC(as.integer(d$X.projid), width = 8, flag = '0', format = 'd')
    d
})

load_phenotypes <- function() load_cohort()[, c(ID_COL, STAGE2_PHENOTYPES)]

# Mean-impute the numeric columns (Figure S2 only; the association engine imputes itself).
mean_impute <- function(df, columns = setdiff(colnames(df), ID_COL)) {
    for (column in columns) {
        if (is.numeric(df[[column]])) df[[column]][is.na(df[[column]])] <- mean(df[[column]], na.rm = TRUE)
    }
    df
}

# ---- Seeded 70/30 split and PGS PCA (Figure 5a-e, S5b, S14) ---------------------------------
# PCA (centred, scaled) on the training rows of the 12 PGS; everyone projected.

fit_pgs_pca <- function() memo('pgs_pca', {
    cohort <- load_cohort()
    pgs_cols <- pgs_variable_ids()
    n <- nrow(cohort)
    set.seed(SPLIT_SEED)
    train_idx <- sample(seq_len(n), size = floor(TRAIN_FRAC * n))
    test_idx <- setdiff(seq_len(n), train_idx)
    pca <- prcomp(cohort[train_idx, pgs_cols], scale. = TRUE)
    list(pca = pca, train_idx = train_idx, test_idx = test_idx,
         variance_explained = summary(pca)$importance[2, ],
         scores_all = predict(pca, cohort[, pgs_cols]), pgs_cols = pgs_cols)
})
