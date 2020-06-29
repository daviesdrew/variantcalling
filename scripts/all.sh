#!/bin/sh

nextflow run ../main.nf \
\
--align all \
\
--variant freebayes --variant_args "-theta freebayes" \
\
--filter bcftools --filter_args "--exclude 'A'" \
\
--reads "../data/CIN*R{1,2}*" \
--ref "../data/ref.fa" \
\
--fastp_min_base_quality 15 \
--fastp_max_percent_low_qual_base 40 \
--verbose 

