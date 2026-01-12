library(ggplot2)
library(glue)
library(ggsci)

plot_pgs_profile <- function(pgs_variables, mapper_pgs, models_directory, plot_directory, chromosome_lengths, mapper_to_rsid, mapper_id_to_csq_group) {
    for (variable in pgs_variables) {
        model_filename <- paste0(models_directory, variable, '.snpnetBETAs.tsv.gz')
        model_df <- read.delim(model_filename, header=TRUE, sep='\t', fill = TRUE, quote="")
        model_df$chrom_pos <- sapply(strsplit(as.character(model_df$X.ID), ":"), function(x) {
            if (x[1] == "XY") {
            x[1] <- "X"
            }
            paste(x[1], x[2], sep=":")
        })
        annotated_path <- paste0(variable, "_tot_unclassified_annon.tsv")
        annotated_supplement <- read.delim(annotated_path, header=TRUE, sep='\t', fill = TRUE, quote="")
        
        model_df$rsid <- sapply(model_df$chrom_pos, mapper_to_rsid)
        model_df$csq_group <- sapply(model_df$chrom_pos, mapper_id_to_csq_group)
        num_unclassified <- sum(model_df$csq_group == "Unclassified")

        print(glue("Number of Unclassified csq_group for {variable}: {num_unclassified}"))
        # Update unclassified variants with information from annotated_supplement
        unclassified_indices <- which(model_df$csq_group == "Unclassified")
        for (idx in unclassified_indices) {
            chrom_pos <- model_df$chrom_pos[idx]
            annotated_match <- annotated_supplement[annotated_supplement$rsid.x == chrom_pos, ]
            if (nrow(annotated_match) > 0) {
                if (!is.na(annotated_match$final_csq[1]) && annotated_match$final_csq[1] != "") {
                    model_df$csq_group[idx] <- annotated_match$final_csq[1]
                }
                if (!is.na(annotated_match$rsid.y[1]) && annotated_match$rsid.y[1] != "") {
                    model_df$rsid[idx] <- annotated_match$rsid.y[1]
                }

            }
        }
        num_na_rsid <- sum(!grepl("^rs", model_df$rsid))
        num_unclassified <- sum(model_df$csq_group == "Unclassified")
        print(glue("Number of Unclassified csq_group for {variable}: {num_unclassified}"))
        write.table(model_df[model_df$csq_group == "Unclassified", ], 
                    file = paste0(plot_directory, gsub("[^A-Za-z0-9_]", "_", mapper_pgs(glue("{variable}_SCORE1_AVG"))), '_unclassified_variants.tsv'), 
                    sep = '\t', row.names = FALSE, quote = FALSE)
        print(glue("Number of unassigned rsid for {variable}: {num_na_rsid}"))
        model_df$chromosome <- sapply(strsplit(as.character(model_df$X.ID), ":"), `[`, 1)
        model_df$position <- as.numeric(sapply(strsplit(as.character(model_df$X.ID), ":"), `[`, 2)) 
        model_df$position <- model_df$position + sapply(model_df$chromosome, chr_to_pos)
        apoe_position_start <- chr_to_pos(19) + 45176340
        apoe_position_end <- chr_to_pos(19) + 45447221
        model_df$top_snps <- FALSE
        top_snps_indices <- order(-abs(model_df$BETA))[1:5]
        model_df$top_snps[top_snps_indices] <- TRUE
        ggplot(model_df, aes(x=position, y=BETA, color=as.factor(sapply(chromosome, chr_to_num) %% 2), shape=csq_group)) +
            geom_point(size = 4) +
            scale_color_manual(values=c("#E69F00", "#56B4E9")) +
            scale_shape_manual(values=c("Intronic"=1, "PAVs"=2, "UTR"=3, "PTVs"=4, "PCVs"=5, "Others"=6, "Unclassified"=7)) +
            theme_bw(base_size=24) +
            labs(
            x = "Genomic position (chromosome)", 
            y = paste0(mapper_pgs(glue("{variable}_SCORE1_AVG")), " PGS weights"),
            shape = "CSQ Group"
            ) +
            geom_hline(yintercept = 0, col = "black") +
            geom_vline(xintercept = c(apoe_position_start, apoe_position_end), linetype="dashed", color="red", size=0.1) +
            scale_x_continuous(
                breaks = sapply(c(1:23), chr_to_pos) + chromosome_lengths[c(1:23)] / 2,
                labels = c(1:9, " ", 11," ", 13, " ", 15," ", 17," ", 19," ", 21," ", "X"),
                limits = c(0, chr_to_pos('MT') + chromosome_lengths[length(chromosome_lengths)])
            ) +
            theme(
            axis.text.x = element_text(angle=90, size=24),
            axis.text.y = element_text(size=24),
            axis.title.x = element_text(size=24),
            axis.title.y = element_text(size=24),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
            ) +
            coord_cartesian(ylim = c(-max(abs(model_df$BETA)), max(abs(model_df$BETA)))) +
            geom_text(
            aes(label = ifelse(model_df$top_snps, model_df$rsid, "")), 
            vjust = 1.5, 
            size = 6,
            )
        ggsave(paste0(plot_directory, gsub("[^A-Za-z0-9_]", "_", mapper_pgs(glue("{variable}_SCORE1_AVG"))), '.pdf'), width = 25, height = 7)
    }


}

chr_to_num <- function(chr) {
    if (chr == "X" || chr == "XY") {
        chr <- 23
    }
    else if (chr == 'Y') {
        chr <- 24
    }
    else if (chr == "MT") {
        chr <- 25
    }
    return(as.numeric(chr))
}

chr_to_pos <- function(chr){
            chr_num <- chr_to_num(chr)
            if (chr_num == 1) {
                return(0)
            }
            return(sum(chromosome_lengths[1:(chr_num-1)]))}