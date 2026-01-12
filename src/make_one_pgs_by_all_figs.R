library(ggplot2)
library(dendextend)
library(dplyr)
library(ggsci)
library(glue)
library(cowplot)
library(patchwork)

plot_pgs_effects <- function(pgs_names, effects, errors, fdr,
  effects_no_apoe, errors_no_apoe, fdr_no_apoe, effects_only_apoe, errors_only_apoe, fdr_only_apoe,
  output_path) {
  
  # Initialize an empty list to store plots
  plot_list <- list()

  effects_df <- rbind(variable = colnames(effects), effects)
  colnames(effects_df) <- NULL
  errors_df <- rbind(variable = colnames(errors), errors)
  colnames(errors_df) <- NULL
  fdr_df <- rbind(variable = colnames(fdr), fdr)
  colnames(fdr_df) <- NULL
  effects_no_apoe_df <- rbind(variable = colnames(effects_no_apoe), effects_no_apoe)
  colnames(effects_no_apoe_df) <- NULL
  errors_no_apoe_df <- rbind(variable = colnames(errors_no_apoe), errors_no_apoe)
  colnames(errors_no_apoe_df) <- NULL
  fdr_no_apoe_df <- rbind(variable = colnames(fdr_no_apoe), fdr_no_apoe)
  colnames(fdr_no_apoe_df) <- NULL
  effects_df <- rbind(effects_df, source = "Full PGS")
  effects_no_apoe_df <- rbind(effects_no_apoe_df, source = "PGS without APOE region")
  errors_df <- rbind(errors_df, source = "Full PGS")
  errors_no_apoe_df <- rbind(errors_no_apoe_df, source = "PGS without APOE region")
  fdr_df <- rbind(fdr_df, source = "Full PGS")
  fdr_no_apoe_df <- rbind(fdr_no_apoe_df, source = "PGS without APOE region")


  for (pgs_name in pgs_names) {
    combined_effects <- cbind(effects_df, effects_no_apoe_df)[c(pgs_name, "variable", "source"), ]
    combined_errors <- cbind(errors_df, errors_no_apoe_df)[pgs_name, ]
    combined_fdr <- cbind(fdr_df, fdr_no_apoe_df)[pgs_name, ]

    # Create a matrix for the asterisks
    asterisks_matrix <- matrix("",
      nrow = nrow(combined_fdr),
      ncol = ncol(combined_fdr))
    asterisks_matrix[
      as.numeric(combined_fdr) < 0.05 & as.numeric(combined_fdr) >= 0.01
    ] <- "0.01 <= FDR < 0.05"
    asterisks_matrix[
      as.numeric(combined_fdr) < 0.01 & as.numeric(combined_fdr) >= 0.001
    ] <- "0.001 <= FDR < 0.01"
    asterisks_matrix[
      as.numeric(combined_fdr) < 0.001
    ] <- "FDR < 0.001"
    asterisks_matrix[
      as.numeric(combined_fdr) >= 0.05
    ] <- "0.05 <= FDR"
    alpha_matrix <- matrix("",
      nrow = nrow(combined_fdr),
      ncol = ncol(combined_fdr))
    alpha_matrix[as.numeric(combined_fdr) > 0.1] <- "0.4"
    alpha_matrix[as.numeric(combined_fdr) <= 0.1] <- "1"

    # Extract the effect sizes and errors for the specified PGS
    pgs_effects <- combined_effects
    pgs_errors <- combined_errors
    pgs_effects <- as.data.frame(t(pgs_effects))
    pgs_errors <- as.data.frame(t(pgs_errors))

    # Create a new dataframe that merges pgs_effects and pgs_errors under new column names
    pgs <- cbind(pgs_effects, pgs_errors)
    colnames(pgs) <- c("effects", "variable", "source", "errors")
    pgs_significance <- as.data.frame(t(asterisks_matrix))
    pgs_alpha <- as.data.frame(t(alpha_matrix))
    pgs$significance <- pgs_significance[, 1]
    pgs$alpha <- pgs_alpha[, 1]
    pgs <- pgs %>% arrange(effects)

    shape_palette <- c("FDR < 0.001" = 16, "0.001 <= FDR < 0.01" = 17,
      "0.01 <= FDR < 0.05" = 18, "0.05 <= FDR" = 22)
    order_vars <- c(
      "Global cognition random slope (demographics)",
      "Global cognitive function (19 tests)",
      "Global AD pathology burden",
      "Neuritic plaque burden (5 regions)",
      "Amyloid level (% cortex area, 8 brain regions)",
      "Diffuse plaque burden (5 regions)",
      "Tangle density (IHC, 8 brain regions)",
      "Neurofibrillary tangle burden (5 regions)",
      "TDP-43 pathology - 4 stages"
      )
    pgs$variable <- factor(pgs$variable, levels = order_vars)

    p <- ggplot(pgs, aes(
      y = as.numeric(effects),
      x = variable,
      color = source,
      group = source,
      shape = significance,
      alpha = alpha     # move alpha mapping here for global use
    )) +
      geom_point(
      size = 3,
      position = position_dodge(width = 0.4)
      ) +
      geom_errorbar(
      aes(
        ymin = as.numeric(effects) - 1.96 * as.numeric(errors),
        ymax = as.numeric(effects) + 1.96 * as.numeric(errors),
        alpha = alpha   # apply same alpha mapping to error bars
      ),
      width = 0.4,
      position = position_dodge(width = 0.4)
      ) +
      scale_alpha_manual(values = c("0.4" = 0.4, "1" = 1)) +
      geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
      theme_bw(base_size = 20) +
      theme(
      axis.text.y        = element_text(angle = 0, hjust = 1),
      axis.text.x        = element_text(angle = 90,vjust = 1),
      panel.grid.major.y = element_line(color = "grey", size = 0.5),
      panel.grid.minor.y = element_line(color = "grey", size = 0.25),
      plot.margin        = unit(c(0, 0, 0, 0), "cm"),
      panel.spacing      = unit(-1.5, "lines")
      ) +
      labs(
      y     = glue("{pgs_name}"),
      x     = "Phenotype (ROSMAP)",
      color = "APOE Region Inclusion",
      shape = "Significance",
      alpha = "FDR Threshold"
      ) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_color_manual(values = pal_npg("nrc")(3)) +
      scale_shape_manual(values = shape_palette) +
      ylim(-0.4, 0.4)

    # Remove x-axis labels for all but the last plot
    if (pgs_name != tail(pgs_names, n = 1)) {
      p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }

    # Add the plot to the list
    plot_list[[pgs_name]] <- p
  }



  only_apoe_effects_df <- rbind(variable = colnames(effects_only_apoe), effects_only_apoe)
  only_apoe_errors_df <- rbind(variable = colnames(errors_only_apoe), errors_only_apoe)
  only_apoe_fdr_df <- rbind(variable = colnames(fdr_only_apoe), fdr_only_apoe)
  colnames(only_apoe_effects_df) <- NULL
  colnames(only_apoe_errors_df) <- NULL
  colnames(only_apoe_fdr_df) <- NULL
  only_apoe_effects_df <- rbind(only_apoe_effects_df, source = "APOE region only")
  only_apoe_errors_df <- rbind(only_apoe_errors_df, source = "APOE region only")
  only_apoe_fdr_df <- rbind(only_apoe_fdr_df, source = "APOE region only")
  only_apoe_effects <- only_apoe_effects_df[c('apoe_genotype', "variable", "source"), ]
  only_apoe_errors <- only_apoe_errors_df[c('apoe_genotype'), ]
  only_apoe_fdr <- only_apoe_fdr_df[c("apoe_genotype"), ]

  asterisks_matrix <- matrix("",
    nrow = nrow(as.matrix(only_apoe_fdr)),
    ncol = ncol(as.matrix(only_apoe_fdr)))
  asterisks_matrix[
    as.numeric(only_apoe_fdr) < 0.05 & as.numeric(only_apoe_fdr) >= 0.01
  ] <- "0.01 <= FDR < 0.05"
  asterisks_matrix[
    as.numeric(only_apoe_fdr) < 0.01 & as.numeric(only_apoe_fdr) >= 0.001
  ] <- "0.001 <= FDR < 0.01"
  asterisks_matrix[
    as.numeric(only_apoe_fdr) < 0.001
  ] <- "FDR < 0.001"
  asterisks_matrix[
    as.numeric(only_apoe_fdr) >= 0.05
  ] <- "0.05 <= FDR"
  alpha_matrix <- matrix("",
    nrow = nrow(as.matrix(only_apoe_fdr)),
    ncol = ncol(as.matrix(only_apoe_fdr)))
  alpha_matrix[as.numeric(only_apoe_fdr) > 0.1] <- ".4"
  alpha_matrix[as.numeric(only_apoe_fdr) <= 0.1] <- "1"

  pgs_effects <- as.data.frame(t(only_apoe_effects))
  pgs_errors <- as.data.frame(t(only_apoe_errors))
  pgs <- cbind(pgs_effects, pgs_errors)
  colnames(pgs) <- c("effects", "variable", "source", "errors")
  pgs_significance <- as.data.frame(t(asterisks_matrix))
  pgs_alpha <- as.data.frame(t(alpha_matrix))
  pgs$significance <- pgs_significance[, 1]
  pgs$alpha <- pgs_alpha[, 1]
  pgs <- pgs %>% arrange(effects)

  order_vars <- c(
    "Global cognition random slope (demographics)",
    "Global cognitive function (19 tests)",
    "Global AD pathology burden",
    "Neuritic plaque burden (5 regions)",
    "Amyloid level (% cortex area, 8 brain regions)",
    "Diffuse plaque burden (5 regions)",
    "Tangle density (IHC, 8 brain regions)",
    "Neurofibrillary tangle burden (5 regions)",
    "TDP-43 pathology - 4 stages"
  )
  pgs$variable <- factor(pgs$variable, levels = order_vars)

  shape_palette <- c("FDR < 0.001" = 16, "0.001 <= FDR < 0.01" = 17,
    "0.01 <= FDR < 0.05" = 18, "0.05 <= FDR" = 22)

  p <- ggplot(pgs, aes(
    y = as.numeric(effects),
    x = variable,
    color = source,
    group = source,
    shape = significance,
    alpha = alpha
  )) +
    geom_point(
      size = 3,
      position = position_dodge(width = 0)
    ) +
    geom_errorbar(
      aes(
        ymin = as.numeric(effects) - 1.96 * as.numeric(errors),
        ymax = as.numeric(effects) + 1.96 * as.numeric(errors),
        alpha = alpha
      ),
      width = 0.4,
      position = position_dodge(width = 0)
    ) +
    scale_alpha_manual(values = c("0.4" = 0.4, "1" = 1)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
    theme_bw(base_size = 20) +
    theme(
      axis.text.y        = element_text(angle = 0, hjust = 1),
      axis.text.x        = element_text(angle = 90, vjust = 1),
      panel.grid.major.y = element_line(color = "grey", size = 0.5),
      panel.grid.minor.y = element_line(color = "grey", size = 0.25),
      plot.margin        = unit(c(0, 0, 0, 0), "cm"),
      panel.spacing      = unit(-1.5, "lines")
    ) +
    labs(
      y     = "APOE region only",
      x     = "Phenotype (ROSMAP variable)",
      color = "APOE Region Inclusion",
      shape = "Significance",
      alpha = "FDR Threshold"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = c("APOE region only" = "red")) +
    scale_shape_manual(values = shape_palette) +
    ylim(-0.4, 0.4)

  plot_list[["APOE_region_only"]] <- p

  # Combine all plots into one using patchwork
  combined_plot <- wrap_plots(plot_list, ncol = 1)

  # Save the combined plot as a PDF with increased width and height
  ggsave(output_path, plot = combined_plot, width = 10, height = 4 * length(pgs_names), units = "in", dpi = 300, limitsize = FALSE)
}