#!/bin/bash

REF=/path/to/reference.fasta
BAM=""
FASTQ1=""
FASTQ2=""

usage() {
	echo "Usage : $0 [ -b BAM NAME ] [ -f FASTQ1 ] [ -r FASTQ2 ]" 1>&2
}

while getopts "b:f:r:" opt ; do
	case ${opt} in
	b)
		BAM=$OPTARG
		;;
	f)
		FASTQ1=$OPTARG
		;;
	r)
		FASTQ2=$OPTARG
		;;
	esac
done

bwa mem -M ${REF} ${FASTQ1} ${FASTQ2} | samtools view -bhT ${REF} | samtools sort -o ${BAM} -
samtools index ${BAM}
