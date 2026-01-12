library(glue)
library(hash)
family <- 'gaussian'
svd_option <- TRUE
assoc_directory <- '~/AD_correlation_subtypes/'
phenotype_df <- read.csv(file=glue('{assoc_directory}continuous_phenotypes.tsv'), sep='\t', check.names=FALSE)
prs_df <- read.csv(file=glue('{assoc_directory}apoe_genotype.tsv'), sep='\t', check.names=FALSE)
covariates_df <- read.csv(file=glue('{assoc_directory}covariates.tsv'), sep='\t', check.names=FALSE)
out_prefix <- glue('{assoc_directory}association_results/')
print(colnames(prs_df))
prs_traits <- colnames(prs_df)[-1]
print(colnames(phenotype_df))
phenotype_traits <- colnames(phenotype_df)[-1]
full_rosmap <- merge(merge(prs_df, phenotype_df, by='#projid'), covariates_df, by='#projid')
write.table(full_rosmap, file=paste(out_prefix, 'full_input_apoe.tsv', sep=''), sep='\t')
print(dim(full_rosmap))
matrix_pval <- matrix(nrow=length(prs_traits), ncol=length(phenotype_traits))
matrix_coef <- matrix(nrow=length(prs_traits), ncol=length(phenotype_traits))
matrix_error <- matrix(nrow=length(prs_traits), ncol=length(phenotype_traits))
row_index <- 0
pvals_for_fdr = c()
for(i in 1:ncol(full_rosmap)){
  full_rosmap[is.na(full_rosmap[,i]), i] <- mean(full_rosmap[,i], na.rm = TRUE)
}
for (prs in prs_traits){
    row_index <- row_index + 1
    col_index <- 0
    for (pheno in phenotype_traits){
            col_index <- col_index + 1
            x <- full_rosmap[[prs]]
            x <- (x-mean(x))/sd(x)
            y <- full_rosmap[[pheno]]
            if (family == "binomial"){
                y <- (y == "True") | (y == 1)
            }
            if (svd_option){
                y <- (y-mean(y))/sd(y)
            }
            
            relation <- glm(y ~ x + PC1_AVG + PC2_AVG + PC3_AVG + PC4_AVG + PC5_AVG + PC6_AVG + PC7_AVG + PC8_AVG + PC9_AVG + PC10_AVG + msex + age_death + age_bl + kronos, data=full_rosmap,
                    family = family, na.action=na.omit)
            matrix_coef[row_index, col_index] <- summary(relation)$coefficients[2,1]
            matrix_pval[row_index, col_index]  <- summary(relation)$coefficients[2,4]
            matrix_error[row_index, col_index] <- summary(relation)$coefficients[2,2]
            pvals_for_fdr = c(pvals_for_fdr, summary(relation)$coefficients[2,4])
    }
}
pvals_fdr <- p.adjust(pvals_for_fdr, method='fdr')
dim(pvals_fdr) <- c(col_index, row_index)
pvals_fdr <- t(pvals_fdr)
pvals_fdr_df = as.data.frame(pvals_fdr)
coef_df <- as.data.frame(matrix_coef)
pvals_df <- as.data.frame(matrix_pval)
error_df <- as.data.frame(matrix_error)
rownames(pvals_fdr_df) <- prs_traits
colnames(pvals_fdr_df) <- phenotype_traits
print(dim(coef_df))
rownames(coef_df) <- prs_traits
colnames(coef_df) <- phenotype_traits
rownames(pvals_df) <- prs_traits
colnames(pvals_df) <- phenotype_traits
rownames(error_df) <- prs_traits
colnames(error_df) <- phenotype_traits
write.table(coef_df, file=paste(out_prefix, 'effects.tsv', sep=''), sep='\t')
write.table(pvals_df, file=paste(out_prefix, 'pvals.tsv', sep=''), sep='\t')
write.table(pvals_fdr_df, file=paste(out_prefix, 'fdr.tsv', sep=''), sep='\t')
write.table(error_df, file=paste(out_prefix, 'errors.tsv', sep=''), sep='\t')
