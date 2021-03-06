library(WadingPool)
PDIR <- '/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/mutect_common/output'
setwd(PDIR)

files <- list.files(PDIR, pattern="processed.vcf$")
samples <- unique(gsub("^chr.*?\\.(.*).vcf$", "\\1", files))
chrs <- paste0("chr", c(1:22, "X", "Y"))

#snp_results <- lapply(samples, genSnpZygosity)

s_paths <- file.path("tmp", gsub("[^a-zA-Z0-9]", "", samples), 
                     paste0("[ID].", samples, "_out.tsv"))

