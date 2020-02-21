#!/bin/bash

########################################
#
# annotate.sh
# Annotate variant information using SnpEff
#
#######################################

snpeff=/synbiodata/thomas/tools/snpEff/snpEff.jar
snpsift=/synbiodata/thomas/tools/snpEff/SnpSift.jar
gnomad=/synbiodata/thomas/AlphaLiquid/bundle/resources/reference/gatk_bundle/af-only-gnomad.hg38.vcf.gz
dbsnp=/synbiodata/thomas/AlphaLiquid/bundle/resources/reference/gatk_bundle/dbsnp_146.hg38.vcf.gz
cosmic=/synbiodata/reference/cosmic/v90/CosmicCodingMuts_v90_chradd.vcf.gz

## VCF filename
VCF=$1
PREFIX=`basename -s ".vcf" ${VCF}`
DIRNAME=`dirname ${VCF}`
ANN=${PREFIX}.snpeff.vcf
DBSNP_ANNO=${PREFIX}.dbsnp.vcf
COSMIC_ANNO=${PREFIX}.dbsnp.cosmic.vcf
GNOMAD_ANNO=${PREFIX}.dbsnp.cosmic.gnomad.vcf
#bcftools sort -o ${DIRNAME}/${PREFIX}.sorted.vcf ${VCF}
java -Xmx4G -Djava.io.tmpdir=/tmp -jar ${snpeff} -noInteraction -noMotif hg38 ${VCF} > ${DIRNAME}/${ANN}
java -Xmx4G -Djava.io.tmpdir=/tmp -jar ${snpsift} annotate -v ${dbsnp} ${DIRNAME}/${ANN} > ${DIRNAME}/${DBSNP_ANNO}
java -Xmx4G -Djava.io.tmpdir=/tmp -jar ${snpsift} annotate -v ${cosmic} ${DIRNAME}/${DBSNP_ANNO} > ${DIRNAME}/${COSMIC_ANNO}
java -Xmx4G -Djava.io.tmpdor=/tmp -jar ${snpsift} annotate -name GNOMAD_ -info AF,AC ${gnomad} ${DIRNAME}/${COSMIC_ANNO} > ${DIRNAME}/${GNOMAD_ANNO}
