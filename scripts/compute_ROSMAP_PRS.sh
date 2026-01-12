#!/bin/bash


ml load plink2

genotype=$1
UKB_PRS_directory=$2
condition=$3
out_directory=$4
apoe=1

if [[ $apoe -eq 1 ]]; then
    plink2 \
        --pfile ${genotype} \
        --silent \
        --out ${out_directory}${condition}.snpnet \
        --score \
        ${UKB_PRS_directory}${condition}.snpnetBETAs.tsv.gz \
        1 2 3 header list-variants
else
    plink2 \
        --pfile ${genotype} \
        --silent \
        --out ${out_directory}${condition}.snpnet \
        --window 1000 \
        --exclude-snp 19:45411941:T:C \
        --score \
        ${UKB_PRS_directory}${condition}.snpnetBETAs.tsv.gz \
        1 2 3 header list-variants
fi
