#!/bin/bash

#SBATCH -t 5:00 
#SBATCH --mem=20G 
#SBATCH -J gatk 
#SBATCH -p all 
#SBATCH -c 4 
#SBATCH -N 1 
#SBATCH -o %x-%j.out 

PATH=$PATH:/mnt/work1/users/pughlab/bin/conda/anaconda3/bin

export PATH
chr='all'
TUMOR="NET-2-011_02_T_DNA.processed"
dbsnpDir='/mnt/work1/users/pughlab/references/dbsnp/dbsnp_b146_GRCh37p13/chr_bed'
dbsnpId='common_all_20151104.bed'
PDIR='/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/gatk_common'

while getopts 'ct:d:f:p:' flag; do
  case "${flag}" in
    c) chr=${OPTARG} ;;
    t) TUMOR=${OPTARG} ;;
    d) dbsnpDir=${OPTARG} ;;
    f) dbsnpId=${OPTARG} ;;
    p) PDIR=${OPTARG} ;;
  esac
done

if($verbose); then
  printf "Argument chr: %s\n" "$chr"
  printf "Argument TUMOR: %s\n" "$TUMOR"
  printf "Argument dbsnpDir: %s\n" "$dbsnpDir"
  printf "Argument dbsnpId: %s\n" "$dbsnpId"
  printf "Argument PDIR: %s\n" "$PDIR"
fi

module load igenome-human/hg19
module load gatk/4.0.5.1

time gatk CollectAllelicCounts \
-I ${PDIR}/input/${TUMOR}.bam \
-R ${REF} \
-L ${dbsnpDir}/${chr}.${dbsnpId} \
-O ${PDIR}/output/${TUMOR}.allelicCounts.tsv