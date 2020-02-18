#!/bin/bash

########################################
#
# annotate.sh
# Annotate variant information using SnpEff
#
#######################################

snpeff=/path/to/snpEff.jar
snpsift=/path/to/SnpSift.jar
gnomad=/path/to/gnomad.vcf.gz
dbsnp=/path/to/dbsnp.vcf.gz
cosmic=/path/to/cosmic_coding_muts.vcf.gz

## VCF filename
VCF=$1
PREFIX=`basename -s ".vcf" ${VCF}`
ANN=${PREFIX}.snpeff.vcf
DBSNP_ANNO=${PREFIX}.dnsnp.vcf
COSMIC_ANNO=${PREFIX}.dbsnp.cosmic.vcf
GNOMAD_ANNO=${PREFIX}.dbsnp.cosmic.gnomad.vcf

java -Xmx4G -Djava.io.tmpdir=/tmp -jar ${snpeff} -noInteraction -noMotif -hgvs1LetterAa hg38 ${VCF} > ${ANN}
java -Xmx4G -Djava.io.tmpdir=/tmp -jar ${snpsift} annotate -v ${dbsnp} ${ANN} > ${DBSNP_ANNO}
java -Xmx4G -Djava.io.tmpdir=/tmp -jar ${snpsift} annotate -v ${cosmic} ${DBSNP_ANNO} > ${COSMIC_ANNO}
java -Xmx4G -Djava.io.tmpdor=/tmp -jar ${snpsift} annotate -name GNOMAD_ -info AF,AC ${gnomad} ${COSMIC_ANNO} > ${GNOMAD_ANNO}
