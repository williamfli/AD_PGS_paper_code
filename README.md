# AD PGS paper code

Code for Li et al., *Cross-trait polygenic scores identify genetic correlates of Alzheimer's
disease heterogeneity*. 

The code reproduces main text figures 2-5, supplementary figures S1-S2, supplementary figures S4-S14, tables S3-S7, and tables S10-S11.

Non-reproduced results: Figure 1 and Figure S3 are schematics. The information for table S1 is available at https://www.radc.rush.edu/docs/documentation.htm. The information for table S2 can be found among the supplementary tables of Tanigawa et al, Significant sparse polygenic risk scores across 813 traits in UK Biobank, PLOS Genetics 2022. Tables S8 and S9 were developed using the web-based Genomic Regions Enrichment of Annotations Tool (GREAT).

The remaining figures and tables can be reproduced through the following steps:

1. Edit scripts/paths.R to assign input and output file paths based on the documentation below.
2. Run "Rscript scripts/run_all.R"


## File format documentation

**`output_dir`** - Output results directory

**`cohort_table`** — TSV with expected columns as follows:
- `#projid` (chr): participant id
- 12 PGS columns named `{GBE_ID}_SCORE1_AVG` (num);
  These 12 should correspond to the prioritized PGS.
- Covariates (num): `PC1_AVG` … `PC10_AVG` (ancestry PCs), `msex` (0/1), `age_bl`,
  `age_death`, `kronos` (genotyping batch, 0 = Broad, 1 = TGen). 
- `apoe_genotype` (num): numeric APOE gradient in {-2, -1, -0.5, 0, 1, 2} for ε2/ε2, ε2/ε3, ε2/ε4,
  ε3/ε3, ε3/ε4, ε4/ε4.
- Phenotypes (num): the 36 Stage 2 variables, titled by phenotype ids

**`all_pgs_scores`** — TSV
- `#projid` (chr): participant id
- 713 PGS columns named `{GBE_ID}_SCORE1_AVG` (num): These are all the PGS tested in Stage 1.

**`pgs_scores_no_apoe`** — TSV
- `#projid` (chr): participant ID 
- 12 PGS columns named `{GBE_ID}_SCORE1_AVG` (num): These should correspond to prioritized PGS that have been rescored without APOE-region variants.

**`ancestry_reference`, `ancestry_tgen`, `ancestry_broad`** — 
- Each are plink2 `--score` output (`.sscore`): `#IID`, `ALLELE_CT`, `NAMED_ALLELE_DOSAGE_SUM`, then `PC1_AVG`, `PC2_AVG`, …. 
- `ancestry_tgen` and `ancestry_broad` should be generated according to instructions here using disease cohort genotypes: https://docs.google.com/document/d/1Z8Vurk49RsTyX9YRhcleXKZomZwG84lj2V_YXjUi1LI/edit?usp=sharing
- `ancestry_reference` is based on Byrska-Bishop, Cell 2022. This file should also have a `SuperPop` column.

**`pgs_annotation`** — Table S1 of Tanigawa et al, PLOS Genetics 2022

**`phenotype_annotation`** — TSV
- `variable` (chr): phenotype id as in `cohort_table`
- `name` (chr): phenotype name for display in the figures

**`chromosome_lengths`** — TSV
- `Chromosome` (chr): 1–22, X, Y, MT in genome order
- `Length` (chr): length of each chromosome in bp, thousands commas allowed

**`great_apob`, `great_prostate`** — TSV export of the GREAT web tool

**`glove_vectors`** — txt file downloaded from https://nlp.stanford.edu/projects/glove/

**`pgs_models_dir`** — directory containing PGS model weights
- File names should be formatted as `{GBE_ID}.snpnetBETAs.tsv.gz` 
- Within each file, columns should include `#ID` (chr) for variant position (`CHROM:POS:REF:ALT`), `A1` (chr) for effect allele, `BETA` (num) for effect size

**`variant_annotation`** — gzipped TSV with one row per variant of the plotted PGS models
- `ID` (chr): variant ID in the format `CHROM:POS:REF:ALT`
- `rsID` (chr): dbSNP rsID
- `Csq_group` (chr): variant annotation (PTVs, PAVs, PCVs, Intronic, UTR, Others)

## Requirements

R 4.4 with ggplot2, patchwork, ggrepel, cowplot, ComplexHeatmap, circlize, dendextend,
pROC, treemap (or treemapify), glue.
