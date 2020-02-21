#!/bin/bash

# Follows
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# Assume BAM file is MarkDuplicated
REF=/synbiodata/reference/gatk_bundle/Homo_sapiens_assembly38.fasta
dbSNP=/synbiodata/reference/gatk_bundle/dbsnp_146.hg38.vcf.gz
GNOMAD=/synbiodata/reference/gatk_bundle/af-only-gnomad.hg38.vcf.gz
GATK=/synbiodata/program/gatk-4.1.4.1
# Assumes recalibrated
usage() {
	echo "Usage : $0 [ -b in.bam ] [ -o out.vcf ] [ -t targetd.bed ]" 1>&2
}


BAM=""
VCF=""
while getopts "b:o:t:" opt ; do
	case $opt in 
	b)
		BAM=${OPTARG}
		;;
	o)
		VCF=${OPTARG}
		;;
	t)
		BED=${OPTARG}
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

echo "Start Mutect2 at :", `date`

if [[ $BED == "" ]]
then
	${GATK}/gatk --java-options "-Xmx2g -Djava.io.tmpdir=/tmp" \
	Mutect2 \
	-R ${REF} \
	-I ${BAM} \
	--max-mnp-distance 0 \
	-O ${VCF}
else
	${GATK}/gatk --java-options "-Xmx2g -Djava.io.tmpdir=/tmp" \
	Mutect2 \
	-R ${REF} \
	-I ${BAM} \
	--max-mnp-distance 0 \
	-L ${BED} \
	-O ${VCF}
fi

