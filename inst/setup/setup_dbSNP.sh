module load tabix/0.2.6

vcffile='common_all_20180418.vcf.gz'
mkdir chr_vcf chr_bed

## Use tabix to subset VCF into chromosomes
zcat ${vcffile} | grep "^#" > header.txt 
for i in $(seq 1 1 22) X Y; do
  tabix ${vcffile} ${i} > tmp
  cat header.txt tmp > chr_vcf/chr${i}.vcf
done
rm tmp header.txt


## Simple perl script to reformat the VCF into a BED file for SNPs that occupy 1 bp
chrcol=0
poscol=1
rsidcol=2
refcol=3
infocol=7

rm chr_bed/all.${vcffile/.vcf.gz/}.bed
time for i in $(seq 1 1 22) X Y; do
  grep -v "^#" chr_vcf/chr${i}.vcf |\
  perl -ne 'my @lines = split("\t", $_);
    my $ref=$lines['${refcol}'];
    if(length($ref)==1){
      my $chr=$lines['${chrcol}'];
      my $end=$lines['${poscol}'];
      my $rs=$lines['${rsidcol}'];
      my @caf=$lines['${infocol}'] =~ m/CAF=(.*?),.*?;/;
      
      print $chr,"\t",
            $end-1, "\t",
            $end, "\t",
            $rs, "\t";
      printf("%.3f\t.\t.\t.\n", 1-$caf[0]);
    }' \
  > chr_bed/chr${i}.bed >> chr_bed/all.${vcffile/.vcf.gz/}.bed
done  
