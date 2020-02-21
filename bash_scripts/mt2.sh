#!/bin/bash

# Follows
# https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# Assume BAM file is MarkDuplicated
REF=
dbSNP=
GNOMAD=
GATK=
# Assumes recalibrated
usage() {
	echo "Usage : $0 [ -p pon.vcf ] [ -b in.bam ] [ -o out.vcf ] [ -t targets.bed ]" 1>&2
}


BAM=""
VCF=""
BED=""
PON=""
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
	p)
		PON=${OPTARG}
		;;
	esac
done

if [[ $BAM == "" ]]
then
	usage
	echo "Input BAM file is required"
	exit 1
fi
if [[ $VCF == "" ]]
then
	usage
	echo "Output VCF is required"
	exit 1
fi

if [[ ${PON} == "" ]]
then
	usage
	echo "Panel Of Normals (PON) is required"
	exit 1
fi

echo "Start Mutect2 at :", `date`

if [[ $BED == "" ]]
then
	${GATK}/gatk --java-options "-Xmx16g -Djava.io.tmpdir=/tmp" \
	Mutect2 \
	-R ${REF} \
	-I ${BAM} \
	-pon ${PON} \
	--dont-use-soft-clipped-bases \
	--germline-resource ${GNOMAD} \
	--f1r2-tar-gz ${VCF}.f1r2.tar.gz \
	-O ${VCF}.tmp.vcf
else
	${GATK}/gatk --java-options "-Xmx16g -Djava.io.tmpdir=/tmp" \
	Mutect2 \
	-R ${REF} \
	-I ${BAM} \
	-O ${VCF}.tmp.vcf \
	-p ${PON} \
	--dont-use-soft-clipped-bases \
	--germline-resource ${GNOMAD} \
	--f1r2-tar-gz ${VCF}.f1r2.tar.gz \
	-L ${BED}

	${GATK}/gatk --java-options "-Xmx16G -Djava.io.tmpdir=/tmp" \
	LearnReadOrientationModel -I ${VCF}.f1r2.targ.gz -O ${VCF}.read.orientation.model.tar.gz

	${GATK}/gatk --java-options "-Xmx16G -Djava.io.tmpdir=/tmp" \
	GetPileupSummaries \
	-I ${BAM} \
	-V ${GNOMAD} \
	-L ${BED} \
fi
