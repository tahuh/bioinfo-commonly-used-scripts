#!/bin/bash

# Follows
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# Assume BAM file is MarkDuplicated
REF=/synbiodata/reference/gatk_bundle/Homo_sapiens_assembly38.fasta
dbSNP=/synbiodata/reference/gatk_bundle/dbsnp_146.hg38.vcf.gz
GNOMAD=/synbiodata/reference/gatk_bundle/af-only-gnomad.hg38.vcf.gz
GATK=/synbiodata/program/gatk-4.1.4.1
# Assumes recalibrated
BAM=$1
VCF=$2
PREFIX=`basename -s ".bam" ${BAM}`
VCF_PREFIX=`basename -s ".vcf.gz" ${VCF}`
DIRNAME=`dirname ${BAM}`
VCF_DIR=`dirname ${VCF}`

