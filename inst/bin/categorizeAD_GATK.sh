#!/bin/bash
ref_col=2
alt_col=3

print_usage() {
  printf "Usage: ..."
}

while getopts 'r:a:o:t:hV' flag; do
  case "${flag}" in
    r) ref_col=${OPTARG} ;;
    a) alt_col=${OPTARG} ;;
    t) tsv="${OPTARG}" ;;
    o) out="${OPTARG}" ;;
    V) verbose=true ;;
    h) print_usage
       exit 1 ;;
  esac
done

if($verbose); then
  printf "Argument vcf: %s\n" "$tsv"
  printf "Argument out: %s\n" "$out"
  printf "Argument ref_col: %d\n" "$ref_col"
  printf "Argument alt_col: %d\n" "$alt_col"
fi


# PDIR='/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/mutect_common/output/tmp/test'
# cp ${PDIR}/chr21.NET-2-001a_02_MT_DNA.processed.allelicCounts.tsv ${PDIR}/x.vcf
# tsv='/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/mutect_common/output/tmp/test/x.vcf'
# file=$1
# fout=${vcf/.vcf/_test.txt}
# format_col=8
# sample_col=10

## Isolate the AD data and categorize into No coverage (0), Ref.Homozygous (1),
# Heterozygous (2), or Alt.Homozygous (3) classifications
start_time=`date +%s`
grep -v "^@" ${tsv} | tail -n+2 | head | \
perl -ne 'use List::Util qw/sum/;
          chomp $_;
          
          ## Isolate the samples Ref and Alt allelic depth
          my @spl=split(/\t/, $_);
          my $ref = $spl['${ref_col}'];
          my $alt = $spl['${alt_col}'];
          
          ## Categorize the depth per SNP
          my $sums = $ref + $alt;
          if($sums == 0){
            print "0\n"; ## NO COVERAGE
          } else {
            my $frac = $ref / $sums;
            if ($frac == 1) {
              print "1\n"; ## REF.HOMOZYGOUS
            } elsif ($frac == 0) {
              print "3\n"; ## ALT.HOMOZYGOUS
            } else {
              print "2\n"; ## HETEROZYGOUS
            }
          }' > \
${out}
end_time=`date +%s`
printf "Runtime: %d\n" `expr $end_time - $start_time`