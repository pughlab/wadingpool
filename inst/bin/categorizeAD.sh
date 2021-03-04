#!/bin/bash
format_col=8
sample_col=10

print_usage() {
  printf "Usage: ..."
}

while getopts 'f:s:o:v:h' flag; do
  case "${flag}" in
    f) format_col=${OPTARG} ;;
    s) sample_col=${OPTARG} ;;
    v) vcf="${OPTARG}" ;;
    o) out="${OPTARG}" ;;
    h) print_usage
       exit 1 ;;
  esac
done

printf "Argument vcf is %s\n" "$vcf"
printf "Argument out is %s\n" "$out"
printf "Argument format_col is %d\n" "$format_col"
printf "Argument sample_col is %d\n" "$sample_col"


# PDIR='/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/mutect_common/output'
# cp ${PDIR}/chr1.NET-2-001a_02_MT_DNA.processed.vcf ${PDIR}/x.vcf
#vcf='/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/mutect_common/output/x.vcf'
# file=$1
# fout=${vcf/.vcf/_test.txt}
# format_col=8
# sample_col=10

echo "Reading VCF from [${vcf}]..."

## From the 'GT:AD:BQ:DP:FA' format, identify the index of 'AD'
AD_IDX=$(grep -v "^#" ${vcf} | head -1 | \
         perl -ne 'my @spl=split(/\t/, $_);
                   my @format_legend = split(/:/, $spl['${format_col}']);
                   my @ad_idx = grep{$format_legend[$_] =~ /^AD$/} 0..$#format_legend;
                   print "@ad_idx\n"')

## Isolate the AD data and categorize into No coverage (0), Ref.Homozygous (1),
# Heterozygous (2), or Alt.Homozygous (3) classifications
time grep -v "^#" ${vcf} | \
perl -ne 'use List::Util qw/sum/;
          chomp $_;
          ## Isolate the samples 'AD' column from 'GT:AD:BQ:DP:FA' data
          my @spl=split(/\t/, $_);
          my @ad_val = split(/:/, $spl['${sample_col}']);
          my @ads = split(/,/, @ad_val['${AD_IDX}']);
          
          ## Categorize the depth per SNP
          my $sums = sum(@ads);
          if($sums == 0){
            print "0\n"; ## NO COVERAGE
          } else {
            my $frac = $ads[0] / $sums;
            if ($frac == 1) {
              print "1\n"; ## REF.HOMOZYGOUS
            } elsif ($frac == 0) {
              print "3\n"; ## ALT.HOMOZYGOUS
            } else {
              print "2\n"; ## HETEROZYGOUS
            }
          }' > \
${out}

echo "... Outputted VCF classifications to [${out}]"