devtools::install_github("quevedor2/WadingPool", ref = "dev")


library(WadingPool)
PDIR <- '/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/mutect_common/output'
setwd(PDIR)

files <- list.files(PDIR, pattern="processed.vcf$")
samples <- unique(gsub("^chr.*?\\.(.*).vcf$", "\\1", files))
chrs <- paste0("chr", c(1:22, "X", "Y"))

snp_results <- lapply(samples, genSnpZygosity)

s_paths <- file.path("tmp", gsub("[^a-zA-Z0-9]", "", samples), 
                     paste0("[ID].", samples, "_out.tsv"))
agg_results <- aggregateAndFilter(s_paths, num_samples=2, dbsnp_file='all.common_all_20151104.bed',
                                  dbsnp_dir=file.path("~", "pughlab/references/dbsnp/dbsnp_b146_GRCh37p13/common/bed")) 

filt_dir <- file.path("tmp", "combined", "filt")
getSampleSimilarity(filt_dir, samples)
