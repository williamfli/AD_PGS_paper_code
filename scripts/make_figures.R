library(ggplot2)
library(reshape2)
library(glue)
library(patchwork)
library(hash)
library(ggrepel)

final_assoc_directory <- "../src/"

source(glue("{final_assoc_directory}make_one_pgs_by_all_figs.R"))
source(glue("{final_assoc_directory}make_main_assoc_figs.R"))
source(glue("{final_assoc_directory}make_pgs_fig.R"))
source(glue("{final_assoc_directory}rGreat_analysis.R"))
source(glue("{final_assoc_directory}transferability_analysis.R"))
source(glue("{final_assoc_directory}make_one_phenotype_by_all.R"))
source(glue("{final_assoc_directory}make_heterogeneity_by_pgs_figs.R"))

phenotype_annotations <- read.delim(glue("{final_assoc_directory}metadata_sample_attributes.tsv"), 
                                    header = TRUE, sep = "\t", fill = TRUE)
pgs_annotations <- read.delim(glue("{final_assoc_directory}journal.pgen.1010105.tsv"), 
                              header = TRUE, sep = "\t", fill = TRUE)
new_header <- as.character(unlist(pgs_annotations[1, ]))
names(pgs_annotations) <- new_header
pgs_annotations <- pgs_annotations[-1, ]
colnames(pgs_annotations)[colnames(pgs_annotations) == "Trait Name"] <- "name"
colnames(pgs_annotations)[colnames(pgs_annotations) == "GBE ID"] <- "variable"
pgs_annotations$variable <- paste0(pgs_annotations$variable, "_SCORE1_AVG")
map_variable_to_name <- function(df) {
  mapping <- setNames(df$name, df$variable)
  map_single_variable <- function(variable) {
    if (variable %in% names(mapping)) {
      return(mapping[[variable]])
    } else {
      return(variable)
    }
  }
  return(map_single_variable)
}
mapper_pheno <- map_variable_to_name(phenotype_annotations)
mapper_pgs <- map_variable_to_name(pgs_annotations)

# Figure for transferability analysis
phenotypes <- read.delim(glue("{final_assoc_directory}continuous_phenotypes.tsv"), header = TRUE, sep = "\t", fill = TRUE, colClasses = c("X.projid" = "character"))
phenotypes <- phenotypes[, colnames(phenotypes) != "tdp_stage4"]
pgs <- read.delim(glue("{final_assoc_directory}combined_prs.tsv"), header = TRUE, sep = "\t", fill = TRUE, colClasses = c("X.projid" = "character"))
shared_ids <- intersect(phenotypes$X.projid, pgs$X.projid)
phenotypes <- phenotypes[phenotypes$X.projid %in% shared_ids, ]
pgs <- pgs[pgs$X.projid %in% shared_ids, ]
full_stage_2_input_for_transferability <- read.delim(glue("{final_assoc_directory}full_input_for_transferability_with_continuous_apoe.tsv"),
                                  header = TRUE, sep = "\t", fill = TRUE, colClasses = c("X.projid" = "character"))
make_transferability_violin(full_stage_2_input_for_transferability, glue("{final_assoc_directory}main_text_figures/transferability/"))

# Figures for stage 1 analysis
stage1_assoc_directory <- glue("{final_assoc_directory}first_round_association_results/")
effects <- read.delim(glue("{stage1_assoc_directory}effects.tsv"), header = TRUE, sep = "\t", fill = TRUE)
errors <- read.delim(glue("{stage1_assoc_directory}errors.tsv"), header = TRUE, sep = "\t", fill = TRUE)
fdr <- read.delim(glue("{stage1_assoc_directory}fdr.tsv"), header = TRUE, sep = "\t", fill = TRUE)
original_pgs_variables <- rownames(effects)
rownames(effects) <- sapply(rownames(effects), mapper_pgs)
colnames(effects) <- sapply(colnames(effects), mapper_pheno)
rownames(errors) <- sapply(rownames(errors), mapper_pgs)
colnames(errors) <- sapply(colnames(errors), mapper_pheno)
rownames(fdr) <- sapply(rownames(fdr), mapper_pgs)
colnames(fdr) <- sapply(colnames(fdr), mapper_pheno)
create_heatmap(effects, fdr,
               glue("{final_assoc_directory}main_text_figures/pgs_effects_heatmap_stage_1_full"))
significant_pgs <- rownames(fdr)[apply(fdr, 1, function(row) any(row < 0.5, na.rm = TRUE))]
effects <- effects[significant_pgs, , drop = FALSE]
errors <- errors[significant_pgs, , drop = FALSE]
fdr <- fdr[significant_pgs, , drop = FALSE]
heatmap_plot_stage_one <- create_heatmap(effects, fdr,
                                         glue("{final_assoc_directory}main_text_figures/pgs_effects_heatmap_stage_1"))

# Main Results starting from stage 2 analysis
apoe_option <- ""
results_directory <- glue("{final_assoc_directory}second_round_association_results_final/")
effects <- read.delim(glue("{results_directory}effects.tsv"), header = TRUE, sep = "\t", fill = TRUE)
errors <- read.delim(glue("{results_directory}errors.tsv"), header = TRUE, sep = "\t", fill = TRUE)
fdr <- read.delim(glue("{results_directory}fdr.tsv"), header = TRUE, sep = "\t", fill = TRUE)
original_pgs_variables <- rownames(effects)
association_selected_phenotypes <- colnames(effects)
rownames(effects) <- sapply(rownames(effects), mapper_pgs)
colnames(effects) <- sapply(colnames(effects), mapper_pheno)
rownames(errors) <- sapply(rownames(errors), mapper_pgs)
colnames(errors) <- sapply(colnames(errors), mapper_pheno)
rownames(fdr) <- sapply(rownames(fdr), mapper_pgs)
colnames(fdr) <- sapply(colnames(fdr), mapper_pheno)
results_directory <- glue("{final_assoc_directory}second_round_association_results_without_apoe_final/")
effects_without_apoe <- read.delim(glue("{results_directory}effects.tsv"), header = TRUE, sep = "\t", fill = TRUE)
errors_without_apoe <- read.delim(glue("{results_directory}errors.tsv"), header = TRUE, sep = "\t", fill = TRUE)
fdr_without_apoe <- read.delim(glue("{results_directory}fdr.tsv"), header = TRUE, sep = "\t", fill = TRUE)
rownames(effects_without_apoe) <- sapply(rownames(effects_without_apoe), mapper_pgs)
colnames(effects_without_apoe) <- sapply(colnames(effects_without_apoe), mapper_pheno)
rownames(errors_without_apoe) <- sapply(rownames(errors_without_apoe), mapper_pgs)
colnames(errors_without_apoe) <- sapply(colnames(errors_without_apoe), mapper_pheno)
rownames(fdr_without_apoe) <- sapply(rownames(fdr_without_apoe), mapper_pgs)
colnames(fdr_without_apoe) <- sapply(colnames(fdr_without_apoe), mapper_pheno)

print(glue("Total number of significant associations {sum(sapply(fdr, function(x) (x < .05)))}"))
print(glue("Dimensions of association matrix: {dim(fdr)}"))

heatmap_plot_stage_two_with_apoe <- create_heatmap(effects, fdr, 
                                                    glue("{final_assoc_directory}main_text_figures/pgs_effects_heatmap_stage_2_final{apoe_option}"))
create_scatter_plot(full_stage_2_input_for_transferability, c("INI30640_SCORE1_AVG", "cancer1044_SCORE1_AVG"), c("plaq_d", "plaq_n"), effects, errors, fdr,
                    mapper_pgs, mapper_pheno,
                    glue("{final_assoc_directory}main_text_figures/scatter_plots/combined_scatter_plots.pdf"))
results_directory <- glue("{final_assoc_directory}second_round_apoe_only/")
effects_just_apoe <- read.delim(glue("{results_directory}effects.tsv"), header = TRUE, sep = "\t", fill = TRUE)
errors_just_apoe <- read.delim(glue("{results_directory}errors.tsv"), header = TRUE, sep = "\t", fill = TRUE)
fdr_just_apoe <- read.delim(glue("{results_directory}fdr.tsv"), header = TRUE, sep = "\t", fill = TRUE)
colnames(effects_just_apoe) <- sapply(colnames(effects_just_apoe), mapper_pheno)
colnames(errors_just_apoe) <- sapply(colnames(errors_just_apoe), mapper_pheno)
colnames(fdr_just_apoe) <- sapply(colnames(fdr_just_apoe), mapper_pheno)

selected_phenotypes <- c(
  "Global cognitive function (19 tests)",
  "Global AD pathology burden",
  "Neuritic plaque burden (5 regions)",
  "Amyloid level (% cortex area, 8 brain regions)",
  "Diffuse plaque burden (5 regions)",
  "Tangle density (IHC, 8 brain regions)"
)
selected_pgs <- c(
  "Alzheimer's/dementia (FH)",
  "Apolipoprotein B",
  "Use of sun/uv protection (Always)",
  "Doctor diagnosed hayfever or allergic rhinitis",
  "Prostate cancer"
)
selected_effects <- effects[, selected_phenotypes, drop = FALSE]
selected_errors <- errors[, selected_phenotypes, drop = FALSE]
selected_fdr <- fdr[, selected_phenotypes, drop = FALSE]
selected_effects_without_apoe <- effects_without_apoe[, selected_phenotypes, 
                                                       drop = FALSE]
selected_errors_without_apoe <- errors_without_apoe[, selected_phenotypes, 
                                                     drop = FALSE]
selected_fdr_without_apoe <- fdr_without_apoe[, selected_phenotypes, 
                                               drop = FALSE]
selected_effects_just_apoe <- effects_just_apoe[, selected_phenotypes, 
                                                 drop = FALSE]
selected_errors_just_apoe <- errors_just_apoe[, selected_phenotypes,
                                                 drop = FALSE]
selected_fdr_just_apoe <- fdr_just_apoe[, selected_phenotypes,
                                               drop = FALSE]
plot_pgs_effects(
  selected_pgs, selected_effects, selected_errors, selected_fdr,
  selected_effects_without_apoe, selected_errors_without_apoe,
  selected_fdr_without_apoe, selected_effects_just_apoe,
  selected_errors_just_apoe, selected_fdr_just_apoe,
  paste0(
    glue("{final_assoc_directory}main_text_figures/one_PGS_many_phenotypes/"),
    "pgs_effects_boxplot_final.pdf"
  ))

# Plot the PGS profiles
variant_annotation_file <- "~/Downloads/ukb_pvar_annotations.gz"
variant_annotations <- read.delim(variant_annotation_file, header = TRUE, sep = "\t", fill = TRUE)
variant_annotations$CHROM_POS <- sapply(strsplit(variant_annotations$ID, ":"), function(x) paste(x[1], x[2], sep = ":"))
id_to_ukb_mapping <- hash::hash(keys = variant_annotations$CHROM_POS, values = variant_annotations$ID_UKB)
id_to_csq_group_mapping <- hash::hash(keys = variant_annotations$CHROM_POS, values = variant_annotations$Csq_group)

map_id_to_ukb <- function(id) {
    if (has.key(id, id_to_ukb_mapping)) {
        return(id_to_ukb_mapping[[id]])
    } else {
        return(id)
    }
}

map_id_to_csq_group <- function(id) {
    if (has.key(id, id_to_csq_group_mapping)) {
        return(id_to_csq_group_mapping[[id]])
    } else {
        return("Unclassified")
    }
}
pgs_models_directory <- "~/Downloads/models/"
chromosome_lengths <- read.delim(glue("{final_assoc_directory}chromosome_lengths.tsv"),
                                header = TRUE, sep = "\t", fill = TRUE)
chromosome_lengths <- as.numeric(gsub(",", "", chromosome_lengths$Length))
pgs_variables <- c('cancer1044', 'INI30640', 'BIN_FC40002267', 'BIN_FC10002267','BIN_FC20002267', 'BIN_FC30002267','FH1263', 
'BIN22126')
plot_pgs_profile(pgs_variables, mapper_pgs, pgs_models_directory,
                glue("{final_assoc_directory}main_text_figures/pgs_illustrations/"), chromosome_lengths, map_id_to_ukb, map_id_to_csq_group)


# Make GREAT enrichment plots
pgs_models_directory <- glue("~/pgs_models_directory/")
plot_pgs_great_enrichment(c('cancer1044', 'INI30640'),
                                pgs_models_directory,
                                glue("{final_assoc_directory}main_text_figures/great_enrichment/"), mapper_pgs)
make_wordmap(glue("{final_assoc_directory}wordmap_data.tsv"))


# Make heterogeneity by PGS figures
selected_phenotype_variables <- c(
  "amyloid", "nft", "plaq_n", "plaq_d", "tangles", "gpath", "cogn_global_lv"
)
train_test_pca <- compare_xgboost_models(
  full_stage_2_input_for_transferability, original_pgs_variables, selected_phenotype_variables, 
  mapper_pgs, mapper_pheno,
  glue("{final_assoc_directory}main_text_figures/heterogeneity_by_pgs/")
)
create_individuals_heatmap(phenotypes, pgs, mapper_pheno, mapper_pgs, full_stage_2_input_for_transferability,original_pgs_variables, train_test_pca$train_idx,
                           glue("{final_assoc_directory}main_text_figures/individuals_heatmap_final.pdf"))
plot_heterogeneity_pca(
  full_stage_2_input_for_transferability, original_pgs_variables, selected_phenotype_variables,
  mapper_pgs, mapper_pheno,
  glue("{final_assoc_directory}main_text_figures/heterogeneity_by_pgs/"),
  train_test_pca$train_idx, train_test_pca$test_idx
)
plot_heterogeneity_pca_phenotypes(
  full_stage_2_input_for_transferability, association_selected_phenotypes,
  mapper_pgs, mapper_pheno,
  glue("{final_assoc_directory}main_text_figures/heterogeneity_by_pgs/")
)
cluster_individuals_by_pgs(
  full_stage_2_input_for_transferability, original_pgs_variables,
  mapper_pgs,
  glue("{final_assoc_directory}main_text_figures/heterogeneity_by_pgs/"),
  train_test_pca$train_idx, train_test_pca$test_idx
)