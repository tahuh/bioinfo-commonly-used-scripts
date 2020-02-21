#!/bin/bash

GATK=/synbiodata/program/gatk-4.1.4.1

usage() {
	echo "Usage : $0 [ -b in.bam ] [ -B out.bam ] [ -M metrics.txt ]" 1>&2 
}

IBAM=""
OBAM=""
METRICS=""
while getopts "b:B:M:" opt; do
	case $opt in
	b)
		IBAM=${OPTARG}
		;;
	B)
		OBAM=${OPTARG}
		;;
	M)
		METRICS=${OPTARG}
		;;
	esac
done

if [[ ${IBAM} == "" ]]
then
	usage
	exit 1
fi

if [[ ${OBAM} == "" ]]
then
	usage
	exit 1
fi

if [[ ${METRICS} == "" ]]
then
	usage
	exit 1
fi

echo "Start Markduplicates at : ", `date`

# ${GATK}/gatk --java-options "-Xmx2G -Djava.io.tmpdir=/tmp" \
# MarkDuplicates \
# I=${IBAM} \
# O=${OBAM} \
# M=${METRICS}

# samtools index ${OBAM}
