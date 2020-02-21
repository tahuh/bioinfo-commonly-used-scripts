#!/bin/bash

REF=
dbSNP=
GNOMAD=
GATK=

BAM=""
OBAM=""
usage() {
	echo "Usage : $0 [ -b in.bam ] [ -o out.bam ] " 1>&2
}

while getopts "b:o:" opt ; do
	case $opt in
	b)
		BAM=${OPTARG}
		;;
	o)
		OBAM=${OPTARG}
		;;
	esac
done

if [[ ${BAM} == "" ]]
then
	usage
	exit 1
fi

if [[ ${OBAM} == "" ]]
then
	usage
	exit 1
fi

BAM_PREFIX=`basename -s ".bam" ${BAM}`
BAM_DIR=`dirname ${BAM}`


echo "Start BQSR at: ",`date`

${GATK}/gatk --java-options "-Xmx2G -Djava.io.tmpdir=/tmp" \
BaseRecalibrator \
-I ${BAM} \
-R ${REF} \
--known-sites ${dbSNP} \
-O ${BAM_DIR}/${BAM_PREFIX}.bqsr.table

${GATK}/gatk --java-options "-Xmx4G -Djava.io.tmpdir=/tmp" \
ApplyBQSR \
-R ${REF} \
-I ${BAM} \
-O ${OBAM} \
--bqsr-recal-file ${BAM_DIR}/${BAM_PREFIX}.bqsr.table

samtools index ${OBAM}
