#!/bin/bash

REF=/synbiodata/thomas/AlphaLiquid/bundle/resources/reference/gatk_bundle/Homo_sapiens_assembly38.fasta
BED=/synbiodata/thomas/AlphaLiquid/bundle/resources/reference/gatk_bundle/Homo_sapiens_assembly38.fasta

usage() {
	echo "Usage : $0 [ -b in.bam ] [ -B target.bed ] [ -v out.vcf ]" 1>&2
}

BAM=""
VCF=""
BED=""

while getopts "b:t:v:" opt ; do
	case $opt in
	b) 
		BAM=${OPTARG}
		;;
	v)
		VCF=${OPTARG}
		;;
	t)
		BED=${OPTARG}
		;;
	esac
done

if [[ ${BAM} == "" ]]
then
	usage
	exit 1
fi

if [[ ${VCF} == "" ]]
then
	usage
	exit 1
fi

if [[ ${BED} == "" ]]
then
	freebayes -f ${REF} -b ${BAM} > ${VCF}
else
	freebayes -f ${REF} -t ${BED} -b ${BAM} > ${VCF}
fi

bgzip ${VCF}
tabix ${VCF}.gz
