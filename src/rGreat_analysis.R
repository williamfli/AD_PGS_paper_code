#library(rGREAT)
library(glue)
library(dplyr)
library(spacyr)       
library(treemap)

plot_pgs_great_enrichment <- function(pgs_variables, models_directory, plot_directory, mapper_pgs) {
    for (pgs in pgs_variables){
        file <- glue("{plot_directory}{pgs}_great_results.tsv")
        tb <- read.table(file, header = TRUE, sep = "\t")
        tb <- tb[tb$HyperFdrQ < 0.05, ]
        tb$fold_change <- as.numeric(tb$RegionFoldEnrich)
        tb$p_value <- as.numeric(tb$BinomFdrQ)
        tb$delabel <- ifelse(tb$Desc %in% head(tb$Desc, 10), tb$Desc, NA)
        print(nrow(tb))
        tb$diffexpressed <- "NO"
        tb$diffexpressed[tb$fold_change > 1.5 & tb$p_value < 0.05] <- "UP"
        print(max(tb$fold_change, na.rm=TRUE))
        print(max(-log10(tb$p_value), na.rm=TRUE)/2)
        p <- ggplot(data = tb, aes(x = fold_change, y = -log10(p_value), col=diffexpressed, label = delabel)) +
            geom_vline(xintercept = c(1.5), col = "black", linetype = 'dashed') +
            geom_hline(yintercept = -log10(0.05), col = "black", linetype = 'dashed') + 
            theme_bw(base_size = 20) +
            theme(
            text = element_text(size = 20, color = "black"),
            axis.text = element_text(size = 20, color = "black"),
            axis.title = element_text(size = 20, color = "black"),
            legend.text = element_text(size = 20, color = "black"),
            legend.title = element_text(size = 20, color = "black"),
            plot.title = element_text(size = 20, color = "black", hjust = 0),
            plot.margin = margin(5.5, 5.5, 5.5, 5.5)
            ) +
            geom_point(size = 2) + 
            scale_color_manual(values = c("grey", "blue"),  
             labels = c("Not significant", "Enriched")) +
            labs(color = 'Severe', 
            x = "Fold-change", y = expression("-log"[10]*"q")) + 
            geom_text_repel(size = 5, color = "black") + 
            annotate("text", x = max(tb$fold_change,na.rm=TRUE), y = -log10(0.05), label = "FDR = 0.05", color = "black", size = 6, hjust = 1) +
            annotate("text", x = 1.5, y = max(-log10(tb$p_value),na.rm=TRUE), label = "1.5 fold-change", angle = 90, color = "black", size = 6, hjust = 0) +
            scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(1, NA), breaks = c(1, pretty(tb$fold_change[tb$fold_change > 1], n = 5)))
        pgs_name <- mapper_pgs(glue("{pgs}_SCORE1_AVG"))
        pgs_name <- gsub("[^[:alnum:]_]", "", pgs_name)
        ggsave(glue("{plot_directory}{pgs_name}_great_enrichment.pdf"), plot = p, width = 12, height = 8, dpi = 300)
    }
}

make_wordmap <- function(wordmap_data, n_clusters = 8) {
  # read input
  df <- read.csv(wordmap_data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  spacy_initialize(model = "en_core_web_sm")
  parsed <- spacy_parse(df$Desc, dependency = FALSE, lemma = FALSE, entity = FALSE)
  tokens <- unique(parsed$token)
  # Get word embeddings using spaCy
  # Read only the vectors for words that exist in tokens
  embeddings <- data.frame()
  connection <- file("~/Downloads/spacy_embeddings.txt", "r")
  while(length(line <- readLines(connection, n = 1)) > 0) {
    parts <- strsplit(line, " ")[[1]]
    word <- parts[1]
    if(word %in% tokens) {
      vector_values <- as.numeric(parts[-1])
      embeddings <- rbind(embeddings, c(word, vector_values))
    }
  }
  close(connection)
  
  # Set proper column names and convert to numeric (except first column)
  if(nrow(embeddings) > 0) {
    colnames(embeddings) <- c("word", paste0("V", 1:(ncol(embeddings)-1)))
    embeddings[,-1] <- lapply(embeddings[,-1], as.numeric)
  }

  # Create phrase vectors by averaging word embeddings
  phrase_vectors <- data.frame()

  for(i in 1:nrow(df)) {
    phrase <- df$Desc[i]
    phrase_parsed <- spacy_parse(phrase, dependency = FALSE, lemma = FALSE, entity = FALSE)
    phrase_tokens <- phrase_parsed$token
    
    # Get embeddings for tokens in this phrase
    phrase_embeddings <- embeddings[embeddings$word %in% phrase_tokens, -1]
    # Filter out common stopwords and punctuation
    phrase_embeddings <- embeddings[embeddings$word %in% phrase_tokens & 
                     !embeddings$word %in% c("of", "to", "-"), -1]
    if(nrow(phrase_embeddings) > 0) {
      # Average the embeddings
      avg_embedding <- colMeans(phrase_embeddings, na.rm = TRUE)
      phrase_vectors <- rbind(phrase_vectors, c(phrase, avg_embedding))
    }
  }

  # Set column names
  if(nrow(phrase_vectors) > 0) {
    colnames(phrase_vectors) <- c("phrase", paste0("V", 1:(ncol(phrase_vectors)-1)))
    phrase_vectors[,-1] <- lapply(phrase_vectors[,-1], as.numeric)
  }

  # Perform clustering on phrase vectors
  phrase_matrix <- as.matrix(phrase_vectors[,-1])
  rownames(phrase_matrix) <- phrase_vectors$phrase
  # Save the phrase matrix
  write.csv(phrase_matrix, file = "phrase_matrix.csv", row.names = TRUE)
  # Read in the phrase matrix
  phrase_matrix <- read.csv("phrase_matrix.csv", row.names = 1)
  phrase_matrix <- as.matrix(phrase_matrix)
  # K-means clustering
  set.seed(123)
  kmeans_result <- kmeans(phrase_matrix, centers = n_clusters, nstart = 25)
  clusters <- kmeans_result$cluster

  # Create data frame for treemap
  treemap_data <- data.frame(
    phrase = rownames(phrase_matrix),
    cluster = as.factor(clusters),
    size = 1  # Equal size for all terms, or you could use frequency/importance
  )

  # Add cluster labels
  cluster_labels <- paste("Cluster", 1:n_clusters)
  treemap_data$cluster_label <- cluster_labels[treemap_data$cluster]

  # Create treemap
  library(RColorBrewer)
  # Save the treemap to PDF
  pdf("/main_text_figures/great_enrichment/treemap_plot.pdf", width = 12, height = 8)
  treemap(treemap_data,
          index = c("cluster_label", "phrase"),
          vSize = "size",
          type = "categorical",
          vColor = "cluster",
          palette = brewer.pal(min(n_clusters, 11), "Set3"),
          title = "Clustered Gene Ontology Terms",
          fontsize.labels = c(0, 16),
          fontcolor.labels = c("transparent", "black"),
          fontface.labels = c(2, 1),
          align.labels = list(c("center", "center"), c("left", "top")),
          overlap.labels = 0.5,
          inflate.labels = FALSE)
  dev.off()


  # Create barplot of top 25 GO terms by -log10(BinomFdrQ) colored by cluster
  # Read df once and process efficiently
  df <- read.csv(wordmap_data, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Merge only necessary columns and immediately filter top 25
  df_subset <- df[, c("Desc", "BinomFdrQ")]
  df_subset$neg_log_q <- -log10(df_subset$BinomFdrQ)
  df_subset <- df_subset[order(df_subset$neg_log_q, decreasing = TRUE), ]
  top_25_indices <- head(order(df_subset$neg_log_q, decreasing = TRUE), 25)
  top_25_df <- df_subset[top_25_indices, ]
  
  # Merge only with top 25 terms
  treemap_subset <- treemap_data[treemap_data$phrase %in% top_25_df$Desc, c("phrase", "cluster_label")]
  top_25 <- merge(top_25_df, treemap_subset, by.x = "Desc", by.y = "phrase", all.x = TRUE)

  # Create barplot
  library(ggplot2)
  p_bar <- ggplot(top_25, aes(x = reorder(Desc, neg_log_q), y = neg_log_q, fill = cluster_label)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3", name = "Cluster") +
    labs(title = "Top 25 Gene Ontology Terms by -log10(BinomFdrQ)",
         x = "GO Process",
         y = "-log10(BinomFdrQ)") +
    theme_bw() +
    theme(axis.text.y = element_text(size = 13),
          axis.text.x = element_text(size = 13))

  # Save barplot
  ggsave("main_text_figures/great_enrichment/top25_go_terms_barplot.pdf", plot = p_bar, width = 10, height = 8, dpi = 300)
}
