#!/bin/bash

#file='x.vcf'
file=$1
fout=${file/.vcf/_test.txt}
format_col=8

echo "Reading VCF from [${file}]..."

## From the 'GT:AD:BQ:DP:FA' format, identify the index of 'AD'
AD_IDX=$(grep -v "^#" ${file} | head -1 | \
         perl -ne 'my @spl=split(/\t/, $_);
                   my @format_legend = split(/:/, $spl['${format_col}']);
                   my @ad_idx = grep{$format_legend[$_] =~ /^AD$/} 0..$#format_legend;
                   print "@ad_idx\n"')

## Isolate the AD data and categorize into No coverage (0), Ref.Homozygous (1),
# Heterozygous (2), or Alt.Homozygous (3) classifications
grep -v "^#" ${file} | \
perl -ne 'use List::Util qw/sum/;
          chomp $_;
          ## Isolate the samples 'AD' column from 'GT:AD:BQ:DP:FA' data
          my @spl=split(/\t/, $_);
          my @ad_val = split(/:/, $spl[10]);
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
${fout}

echo "... Outputted VCF classifications to [${fout}]"