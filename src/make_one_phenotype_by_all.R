library(ggplot2)
library(dendextend)
library(dplyr)
library(ggsci)
library(glue)


plot_phenotype_effects <- function(phenotype_name, effects, errors, fdr, output_path, 
                             plot_height, apoe_removed=""){


  # Create a matrix for the asterisks
  asterisks_matrix <- matrix("", nrow = nrow(fdr), ncol = ncol(fdr))
  asterisks_matrix[fdr < 0.05 & fdr >= 0.01] <- "0.01 <= FDR < 0.05"
  asterisks_matrix[fdr < 0.01 & fdr >= 0.001] <- "0.001 <= FDR < 0.01"
  asterisks_matrix[fdr < 0.001] <- "FDR < 0.001"
  asterisks_matrix[fdr >= 0.05] <- "0.05 <= FDR"
  rownames(asterisks_matrix) <- rownames(fdr)
  colnames(asterisks_matrix) <- colnames(fdr)

  # Extract the asterisks for the specified phenotype
  phenotype_asterisks <- asterisks_matrix[, phenotype_name]
  
  # Extract the effect sizes and errors for the specified phenotype
  phenotype_effects <- as.data.frame(effects[, phenotype_name])
  phenotype_errors <- as.data.frame(errors[, phenotype_name])
  
  # Create a new dataframe that merges phenotype_effects and phenotype_errors under new column names
  phenotype <- cbind(phenotype_effects, phenotype_errors)
  colnames(phenotype) <- c("effects", "errors")
  phenotype$variable <- rownames(effects)
  print(phenotype)
  
  # Add the cluster and notes columns to the dataframe
  shape_palette <- c("FDR < 0.001" = 16, "0.001 <= FDR < 0.01" = 17, 
                     "0.01 <= FDR < 0.05" = 18, "0.05 <= FDR" = 22)
p <- ggplot(phenotype, aes(y = reorder(variable, effects), x = effects, 
                                         shape = factor(phenotype_asterisks, levels = c("FDR < 0.001", 
                                           "0.001 <= FDR < 0.01", 
                                           "0.01 <= FDR < 0.05", 
                                            "0.05 <= FDR")))) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin = effects - errors, xmax = effects + errors), width = 0.2) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
    theme_bw(base_size = 16) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1),
                panel.grid.major.y = element_line(color = "grey", size = 0.5),
                panel.grid.minor.y = element_line(color = "grey", size = 0.25)) +
    labs(
             y = "PGS", x = glue("{phenotype_name} Effect size [68% CI] {apoe_removed}"), 
             shape = "Significance") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = pal_npg("nrc")(3)) +
    scale_shape_manual(values = shape_palette) +
    xlim(-0.3, 0.3)

  # Save the plot with increased height and adjust the y-axis text size
  ggsave(output_path, plot = p + theme(axis.text.y = element_text(size = 8)), 
         width = 12, height = plot_height, units = "in", dpi = 300)

  # Save the plot as a PDF
  pdf_output_path <- sub("\\.png$", ".pdf", output_path)
  ggsave(pdf_output_path, plot = p + theme(axis.text.y = element_text(size = 8)), 
         width = 12, height = plot_height, units = "in", dpi = 300)
  return(p)
}