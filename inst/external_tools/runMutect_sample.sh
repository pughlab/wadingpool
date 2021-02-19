#!/bin/bash
#
#$ -cwd

TUMOR="NET-2-001a_02_MT_DNA.processed"
dbsnpDir='/mnt/work1/users/pughlab/references/dbsnp/dbsnp_b146_GRCh37p13/common/bed'
dbsnpId='common_all_20151104.bed'
PDIR='/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/mutect_common'

module load igenome-human/hg19
module load  mutect/1.1.5
module load gatk/3.5
module unload java/6
module load java/7

java -Xmx13g -jar $mutect_dir/muTect.jar \
--analysis_type  MuTect \
--min_qscore 20 \
-R ${REF} \
-I:tumor ${PDIR}/input/chr_subset/chr22/chr22.${TUMOR}.bam \
-L ${dbsnpDir}/chr22.${dbsnpId} \
--force_output \
-vcf ${PDIR}/output/chr22.${TUMOR}.vcf