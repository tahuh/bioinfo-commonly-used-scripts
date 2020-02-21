#!/bin/bash

##############################################
#
# default_pipeline.sh
#
# The basic pipeline of genomics starts from 
# adatpter trimming to analysis read bam
#
# Author : Thomas Sunghoon Heo
#
##############################################

REF=/path/to/indexed/reference/genome
picard=/path/to/picard/jar/file
gatk=/path/to/gatk/jar

## Adapter trimming
fastp \
-i read1.fastq.gz \
-I read2.fastq.gz \
-o read1.trimmed.fastq.gz \
-O read2.trimmed.fastq.gz \
-h sample.html \
-j sample.json \
-q 20 -u 20 -x -y -3 -p -g -t 1 -T 1 \
--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

## Mapping and assign read group id for further analysis
bwa mem -M $REF read1.trimmed.fastq.gz read2.trimmed.fastq.gz | \
samtools view -bhT $REF | \
samtools sort | \
java -Xmx4G -Djava.io.tmpdir=/tmp -jar ${picard} AddOrReplaceReadGroups \
I=/dev/stdin \
O=${sample}.bam \
RGID=${sample} \
RGLB=${sample} \
RGPL=ILLUMINA \
RGPU=SequencingMachine \
RGSM=${sample} \
RGCN=Yonsei \
RGPM=NovaSeq

samtools index ${sample}.bam

java -Xmx4G -Djava.io.tmpdir=/tmp -jar ${picard} MarkDuplicates \
I=${sample}.bam \
O=${sample}.markdup.bam \
M=${sample}.markdup.metrics

samtools index ${sample}.markdup.bam

