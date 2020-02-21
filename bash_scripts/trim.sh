#!/bin/bash

REF=/synbiodata/reference/gatk_bundle/Homo_sapiens_assembly38.fasta
FASTQ1=""
FASTQ2=""
OUT1=""
OUT2=""
usage() {
	echo "Usage : $0 [ -f FASTQ1 ] [ -F FASTQ2 ] [ -o FASTQ1.trimmed ] [ -O FASTQ2.trimmed ]" 1>&2
}

while getopts "f:F:o:O:" opt ; do
	case ${opt} in
	f)
		FASTQ1=$OPTARG
		;;
	F)
		FASTQ2=$OPTARG
		;;
	o)
		OUT1=$OPTARG
		;;
	O)
		OUT2=$OPTARG
		;;
	esac
done

fastp \
-i ${FASTQ1} \
-I ${FASTQ2} \
-o ${OUT1} \
-O ${OUT2} \
-q 20 -u 20 -x -y -3 -p -g -t 1 -T 1 \
--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
