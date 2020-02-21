#!/bin/bash
REF=/synbiodata/reference/gatk_bundle/Homo_sapiens_assembly38.fasta
dbSNP=/synbiodata/reference/gatk_bundle/dbsnp_146.hg38.vcf.gz
GNOMAD=/synbiodata/reference/gatk_bundle/af-only-gnomad.hg38.vcf.gz
GATK=/synbiodata/program/gatk-4.1.4.1
PONDIR=/synbiodata/thomas/Clinical/PON/200207_PON_WITH_PBMCS
PONVCFDIR=${PONDIR}/vcfs
BED=/synbiodata/thomas/AlphaLiquid/bundle/resources/bed/hg38/al100_probe_hg38.bed
echo "Start prepare Panel Of Normals(PON) at: ,",`date`

${GATK}/gatk --java-options "-Xmx32G -Djava.io.tmpdir=/tmp -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
GenomicsDBImport \
-R ${REF} \
-L ${BED} \
--genomicsdb-workspace-path ${PONDIR} \
-V ${PONVCFDIR}/Br_1_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_2_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_3_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_4_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_5_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_6_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_7_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_8_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_9_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_10_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_11_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_12_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_13_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_17_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_18_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_19_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_20_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_21_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_22_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_24_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_25_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_26_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_27_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_28_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_31_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_32_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_33_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_36_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_37_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_38_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_39_pbmc.mt2.vcf \
-V ${PONVCFDIR}/Br_40_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_13_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_14_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_15_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_16_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_17_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_18_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_19_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_20_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_21_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_22_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_23_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_24_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_25_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_26_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_27_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_28_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_29_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_30_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_31_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_32_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_33_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_34_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_36_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_37_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_38_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_39_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_43_pbmc.mt2.vcf \
-V ${PONVCFDIR}/OE_45_pbmc.mt2.vcf

${GATK}/gatk \
CreateSomaticPanelOfNormals -R ${REF} --germline-resource ${GNOMAD} -V genedb://${PONDIR} -O ${PONDIR}/OE_BR_200207_NovaSeq_PBMC_PON.vcf
