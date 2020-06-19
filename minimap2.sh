#!/bin/sh

nextflow run main.nf \
\
--align minimap2 --align_args "-minr_kmer_len minimap2" \
\
--variant freebayes --variant_args "-theta freebayes" \
\
--filter bcftools --filter_args "-exclude 'A'" \
\
--prediction snpeff --prediction_args "-soft_filter snpeff" \
\
--consensus bcftools --consensus_args "-key bcftools" \
\
--reads "./data/CIN*R{1,2}*" \
--ref "./data/ref.fa" \
\
--fastp_min_base_quality 15 \
--fastp_max_percent_low_qual_base 40 \
--verbose 

