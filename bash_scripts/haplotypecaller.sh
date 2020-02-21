#!/bin/bash

# Follows
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# Assume BAM file is MarkDuplicated
REF=/synbiodata/reference/gatk_bundle/Homo_sapiens_assembly38.fasta
dbSNP=/synbiodata/reference/gatk_bundle/dbsnp_146.hg38.vcf.gz
GNOMAD=/synbiodata/reference/gatk_bundle/af-only-gnomad.hg38.vcf.gz
GATK=/synbiodata/program/gatk-4.1.4.1
BAM=""
VCF=""
usage() {
	echo "Usage: $0 [ -b in.bam ] [ -v out.vcf.gz ]" 1>&2
}

while getopts "b:v:" opt ; do
	case $opt in
	b)
		BAM=${OPTARG}
		;;
	v)
		VCF=${OPTARG}
		;;
	esac
done

if [[ $BAM == "" ]]
then
	usage
	exit 1
fi

if [[ $VCF == "" ]]
then
	usage
	exit 1
fi

VCF_PREFIX=`basename -s ".vcf.gz" ${VCF}`
VCF_DIR=`dirname ${VCF}`

echo "Start HaplotypeCaller at: ",`date`

${GATK}/gatk --java-options "-Xmx2g -Djava.io.tmpdir=/tmp" \
HaplotypeCaller \
-R ${REF} \
-I ${BAM} \
-O ${VCF} \
-ERC GVCF


### Below is not working correctly
# echo "Start CNNScoreVariants at: ",`date`

# ${GATK}/gatk --java-options "-Xmx2g -Djava.io.tmpdir=/tmp" \
# CNNScoreVariants \
# -R ${REF} \
# -V ${VCF} \
# -I ${BAM} \
# -O ${VCF_DIR}/${VCF_PREFIX}.cnn.vcf.gz \
# -tensor-type read_tensor

# echo "Start FilterVariantTranches at: ",`date`

# ${GATK}/gatk --java-options "-Xmx2g -Djava.io.tmpdir=/tmp" \
# FilterVariantTranches \
# -V ${VCF_DIR}/${VCF_PREFIX}.cnn.vcf.gz \
# --resource ${dbSNP} \
# --resource ${GNOMAD} \
# --info-key CNN_2D \
# --snp-tranche 99.95 \
# --indel-tranche 99.4 \
# --invalidate-previous-filters \
# -O ${VCF_DIR}/${VCF_PREFIX}.cnn.filter.variant.tranches.vcf.gz
