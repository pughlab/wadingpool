library(GenomicRanges)
library(VariantAnnotation)
library(optparse)

option_list = list(
  make_option(c("-r", "--refdir"), type="character", default=NULL, 
              help="Directory containing the dbSNP common_all file", metavar="character"),
  make_option(c("-v", "--vcf"), type="character", default='common_all_20151104.vcf.gz', 
              help="VCF file for dbSNP common_all", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# opt$refdir <- "/mnt/work1/users/pughlab/references/dbsnp/dbsnp_b151_GRCh38p7/chr_vcf"
# opt$vcf <- 'chr1.vcf'
setwd(opt$refdir)
dir.create(file.path("common", "bed"), recursive = T, showWarnings = F)

standardizeSeqLevels <- function(gr){
  # Assign the UCSC chromosome labelling
  ucsc_style <- mapSeqlevels(seqlevels(gr), "UCSC") 
  gr <- renameSeqlevels(gr, ucsc_style[!is.na(ucsc_style)])
  seqlevels(gr) <- seqlevelsInUse(gr)
  gr
}

dbsnp_common <- readVcf(opt$vcf)
dbsnp_common <- standardizeSeqLevels(dbsnp_common)
for(each_chr in levels(seqnames(dbsnp_common))){
  print(paste0("Formatting ", each_chr, "..."))
  dbsnp_chr <- dbsnp_common[seqnames(dbsnp_common) == each_chr]
  seqlevels(dbsnp_chr) <- seqlevelsInUse(dbsnp_chr)
  
  df <- data.frame(seqnames=seqnames(dbsnp_chr),
                   starts=as.integer(start(dbsnp_chr)-1),
                   ends=as.integer(end(dbsnp_chr)),
                   rs=gsub("^", "rs", info(dbsnp_chr)$RS),
                   caf=sapply(info(dbsnp_chr)$CAF, function(x) 1 - as.numeric(x[1])),
                   strands=strand(dbsnp_chr),
                   names=c(rep(".", length(dbsnp_chr))),
                   scores=c(rep(".", length(dbsnp_chr))))
  write.table(df, file=file.path("common", "bed", 
                                 paste0(each_chr, ".common_all_20151104.bed")), 
              quote=F, sep="\t", row.names=F, col.names=F)
}

