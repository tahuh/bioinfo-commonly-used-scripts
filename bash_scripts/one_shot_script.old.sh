#!/bin/bash

### CONSTANT FILES
REF=
BED=
dbSNP=
gnomad=

### PROGRAMS
gatk=/synbiodata/program/gatk-4.1.4.1

usage() {
	echo "Usage: $0 [ -f read1.fastq.gz ] [ -r read2.fastq.gz ] [ -o output path ] [ -s sample ] [ -t threads ]" 1>&2
}

FASTQ1=""
FASTQ2=""
ODIR=""
SAMPLE=""
THREADS=2
while getopts "f:r:o:s:t:" opt; do
	case $opt in
	f)
		FASTQ1=${OPTARG}
		;;
	r)
		FASTQ2=${OPTARG}
		;;
	o)
		ODIR=${OPTARG}
		;;
	s)
		SAMPLE=${OPTARG}
		;;
	t)
		THREADS=${OPTARG}
		;;
	esac
done

if [[ ${FASTQ1} == "" ]]
then
	usage
	exit 1
fi

if [[ ${FASTQ2} == "" ]]
then
	usage
	exit 1
fi

if [[ ${ODIR} == "" ]]
then
	usage
	exit 1
fi

if [[ ${SAMPLE} == "" ]]
then
	usage
	exit 1
fi


echo "./one_shot_script.sh -f ${FASTQ1} -r ${FASTQ2} -o ${ODIR} -s ${SAMPLE} -t ${THREADS}"



echo "Start fastq trimmind at : ", `date`

s=${SAMPLE}
fastp \
-i ${FASTQ1} \
-I ${FASTQ2} \
-o ${ODIR}/trimmed/${s}_R1.fq.gz \
-O ${ODIR}/trimmed/${s}_R2.fq.gz \
-h ${ODIR}/trimmed/${s}.html \
-j ${ODIR}/trimmed/${s}.json \
-q 20 -u 20 -x -y -3 -p -g -t 1 -T 1 \
--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-w ${THREADS}

echo "Start Alignment at : " , `date`
bwa mem -M -t ${THREADS} ${REF} ${ODIR}/trimmed/${s}_R1.fq.gz ${ODIR}/trimmed/${s}_R2.fq.gz | samtools view -bhT ${REF} | samtools sort -o ${ODIR}/bam/${s}.tmp.bam
samtools index ${ODIR}/${bam}/${s}.tmp.bam

sample=s

echo "Start adding read groups t : ", `date`
${gatk}/gatk --java-options "-Xmx16G -Djava.io.tmpdir=${ODIR}/gatk_tmp" \
AddOrReplaceReadGroups \
-I ${ODIR}/bam/${sample}.tmp.bam \
-O ${ODIR}/bam/${sample}.bam \
RGID=${sample} \
RGLB=${sample} \
RGPL=ILLUMINA \
RGPU=DKDL202001359 \
RGSM=${sample} \
RGCN=Novogen \
RGPM=NovaSeq

samtools index ${ODIR}/bam/${sample}.bam

rm ${ODIR}/bam/${sample}.tmp.bam
rm ${ODIR}/bam/${sample}.tmp.bam.bai


echo "Start pileup at : ", `date`

BAM=${ODIR}/bam/${sample}.bam
samtools mpileup -d 100000 -x -q 0 -Q 0 -f ${REF} ${BAM} > ${ODIR}/pileup/${s}.pileup

echo "Start Markduplicates at : ", `date`
${gatk}/gatk --java-options "-Xmx16G -Djava.io.tmpdir=/tmp" \
MarkDuplicates \
-I ${ODIR}/bam/${s}.bam \
-O ${ODIR}/bam/${s}.markdup.bam \
-M ${ODIR}/bam/${s}.markdup.metrics.txt

samtools index ${ODIR}/bam/${s}.markdup.bam

echo "Start BQSR at: ",`date`
${gatk}/gatk --java-options "-Xmx16G -Djava.io.tmpdir=/tmp" \
BaseRecalibrator \
-I ${ODIR}/bam/${s}.markdup.bam \
-R ${REF} \
--known-sites ${dnSNP} \
-O ${ODIR}/bam/${s}.bqsr.table

${gatk}/gatk --java-options "-Xmx16G -Djava.io.tmpdir=/tmp" \
ApplyBQSR \
-R ${REF} \
-I ${ODIR}/bam/${s}.markdup.bam \
-O ${ODIR}/bam/${s}.recal.bam \
--bqsr-recal-file ${ODIR}/bam/${s}.bqsr.table

samtools index ${ODIR}/bam/${s}.recal.bam
