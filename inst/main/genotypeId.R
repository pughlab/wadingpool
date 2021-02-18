# Genotype identifier for shallow WGS
library(GenomicRanges)
library(VariantAnnotation)


pdir <- args[1]
vcfFile <- args[2]

pdir <- '/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/mutect_common'
vcfFile <- 'merge.vcf'

mvcf <- readVcf(file.path(pdir, "output", "filt_all", vcfFile))

genotype <- geno(mvcf)$AD
ignore_idx <- grep("none", colnames(genotype))
ignore_idx <- c(1, 3, 5, 8)
genotype <- genotype[,-ignore_idx]
geno_vals <- apply(genotype, 2, function(sample_vcf){
  zyg.vec <- sapply(sample_vcf, function(x){
    zyg <- NA
    if(length(x) > 0){
      af_id <- try( x / sum(x))
      if(any(is.nan(af_id))) af_id <- c(-1, -1)

      if(af_id[1] == 1) zyg <- 'REF.HOM' else if(af_id[1] > 0) zyg <- 'HET' else if(af_id[1] == -1) zyg <- 'NOCOV' else zyg <- 'ALT.HOM'
    }
    zyg
  })
  zyg.vec
})

colnames(geno_vals) <- c("NET-001a-MT", "NET-001b", "NET-005", "NET-011")
jaccard <- function(gv){
  a_na <- which(is.na(gv[,1]))
  b_na <- which(is.na(gv[,2]))
  
  ab_na <- unique(sort(c(a_na, b_na)))
  gv <- gv[-ab_na,]
  ref_hom_idx <- which(apply(gv, 1, function(x) all(x=='REF.HOM')))
  gv <- gv[-ref_hom_idx,]
  print(paste0("Comparing ", paste(colnames(gv), collapse="-"), "with ", nrow(gv), "rows"))
  vals <- apply(gv, 1, function(x) length(unique(x)))
  length(which(vals == 1)) / length(vals)
}

jac_mat <- sapply(colnames(geno_vals), function(a){
  sapply(colnames(geno_vals), function(b){
    jaccard(geno_vals[,c(a,b)])
  })
})

