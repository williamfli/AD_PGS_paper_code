library(ggplot2)
library(pROC)
library(glue)
library(patchwork)

make_transferability_violin <- function(df,
                                        plot_filename,
                                        width = 10,
                                        height = 10) {

    # NIA Reagan score transferability with FH1263
    df_plot <- df[!is.na(df$niareagansc) & !is.na(df$FH1263_SCORE1_AVG), ]
    print(glue("Total individuals: {nrow(df)}"))
    print(glue("Individuals with complete data: {nrow(df_plot)}"))
    print(glue("Individuals excluded: {nrow(df) - nrow(df_plot)}"))
    df_plot$niareagansc_binary <- ifelse(df_plot$niareagansc %in% c(1, 2), "Case (NIA Reagan Score ≤ 2)", "Control (NIA Reagan Score > 2)")
    df_plot$niareagansc_binary <- as.factor(df_plot$niareagansc_binary)

    auroc <- auc(roc(cases=df_plot$FH1263_SCORE1_AVG[df_plot$niareagansc_binary == "Case (NIA Reagan Score ≤ 2)"], controls=df_plot$FH1263_SCORE1_AVG[df_plot$niareagansc_binary == "Control (NIA Reagan Score > 2)"]))
    print(glue("AUROC niareagansc: {auroc}"))
    ad_median <- median(df_plot$FH1263_SCORE1_AVG[df_plot$niareagansc_binary == "Case (NIA Reagan Score ≤ 2)"])
    not_ad_median <- median(df_plot$FH1263_SCORE1_AVG[df_plot$niareagansc_binary == "Control (NIA Reagan Score > 2)"])
    df_plot$FH1263_SCORE1_AVG <- scale(df_plot$FH1263_SCORE1_AVG)
    
    p1 <- ggplot(df_plot, aes(x=niareagansc_binary, y=FH1263_SCORE1_AVG)) +
        geom_violin(trim = FALSE) +
        geom_boxplot(width = 0.1) +
        theme_bw(base_size = 25) +
        theme(
          axis.text.x        = element_text(angle = 45, hjust = 1, size = 20),
          panel.grid.major.x = element_line(color = "grey", size = 0.5),
          panel.grid.minor.x = element_line(color = "grey", size = 0.25),
          plot.title         = element_text(hjust = 0.5, size = 22)
        ) +
        labs(
          y = "AD Family History PGS (z-score)",
          x = "Alzheimer's disease status based on NIA Reagan Score",
          title = "Violin/Boxplot"
        ) +
        scale_fill_manual(values = rainbow(10))
    
    # compute ROC object
    roc_obj <- roc(
      cases    = df_plot$FH1263_SCORE1_AVG[df_plot$niareagansc_binary == "Case (NIA Reagan Score ≤ 2)"],
      controls = df_plot$FH1263_SCORE1_AVG[df_plot$niareagansc_binary == "Control (NIA Reagan Score > 2)"]
    )
    
    # ROC plot with AUROC annotation
    p2 <- ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, colour = "#1c61b6") +
      theme_bw(base_size = 25) +
      labs(
        title = "ROC curve for NIA Reagan",
        x     = "1 - Specificity",
        y     = "Sensitivity"
      ) +
      annotate(
        "text", x = 0.6, y = 0.2,
        label = glue("AUROC = {round(auroc, 3)}"),
        size  = 6
      ) +
      theme(plot.title = element_text(hjust = 0.5, size = 22))
    
    # save both plots
    ggsave(glue('{plot_filename}niareagansc_transferability.pdf'),
           plot = p1, width = width, height = height)
    ggsave(glue('{plot_filename}niareagansc_ROC.pdf'),
           plot = p2, width = width, height = height)
}