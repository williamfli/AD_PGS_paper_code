# AD_PGS_paper_code

Code for Li et al., *Cross-trait polygenic scores identify genetic correlates of Alzheimer's
disease heterogeneity*. 

The code reproduces main text figures 2-5, supplementary figures S1-S2, supplementary figures S4-S10, tables S3-S7, and tables S10-S12.


Non-reproduced results: Figure 1 and Figure S3 are schematics. The information for table S1 is available at https://www.radc.rush.edu/docs/documentation.htm. The information for table S2 can be found among the supplementary tables of Tanigawa et al, Significant sparse polygenic risk scores across 813 traits in UK Biobank, PLOS Genetics 2022. Tables S8 and S9 were developed using the web-based Genomic Regions Enrichment of Annotations Tool (GREAT).

The remaining tables, and all of the figure panels, can be reproduced through the following steps:

1. Edit scripts/paths.R to assign input and output file paths based on the documentation below.
2. Run "Rscript scripts/run_all.R"


## Documentation for file paths
| File path| Description |
|---|---|
| `output_dir` | where results are written |
| `cohort_table` | individual-level prioritized PGS, AD phenotypes, covariates, and APOE allelotype |
| `all_pgs_scores` | individual-level PGS (all 713) |
| `pgs_scores_no_apoe` | individual-level prioritized PGS, scored without the APOE region |
| `ancestry_reference`, `ancestry_tgen`, `ancestry_broad` | `.sscore` projections of ROSMAP individual genotypes onto PC loadings from the 1000 Genomes Project and Human Genome Diversity Project |
| `pgs_annotation` | Supplementary table from Tanigawa et al, PLOS 2022 (Table S2 in this paper)|
| `phenotype_annotation` | Annotated ROSMAP phenotype names that can be obtained from https://www.radc.rush.edu/docs/documentation.htm. (Table S1)|
| `chromosome_lengths` | hg19 chromosome lengths |
| `great_apob`, `great_prostate` | GREAT enrichment tables for ApoB and prostate cancer PGS obtained using the GREAT web tool (Table S8-S9)| 
| `glove_vectors` | GloVe 2024 Wikipedia+Gigaword 50-d word vectors |
| `pgs_models_dir`, `variant_annotation`| UKB snpnet weight files and variant annotation |


## Requirements

R 4.4 with ggplot2, patchwork, ggrepel, cowplot, ComplexHeatmap, circlize, dendextend,
cluster, pROC, treemap (or treemapify), glue.
