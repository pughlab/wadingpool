# devtools::install_github("quevedor2/WadingPool", ref = "dev")
suppressPackageStartupMessages(library(optparse))
library(WadingPool)

####################
#### Parameters ####
option_list <- list(
  make_option(c("-s", "--sample"), type="character",
              help="Path to the file containing samples genotype data; csv for autosome, vcf for chrM"),
  make_option(c("-i", "--ids"), type="character", default='NA',
              help="Comma-separated list of ordered sample IDs, mainly for chrM"),
  make_option(c("-m", "--mode"), type="character", default='autosome',
              help="Sample similarity for 'autosome' or 'chrM' [default=%default]"),
  make_option(c("-p", "--midpoint"), type="numeric", default=0.5,
              help="Numeric value for midpoint of sample similaritys (Recommended 0.03 for autosome, 0.5 for chrM) [default=%default]"),
  make_option(c("-o", "--outdir"), type="character", default='results/sampleid',
              help="Path to out directory [default=%default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

sample    <- opt$sample     # 
mode      <- opt$mode       # 
ids       <- opt$ids       # 
midpoint  <- opt$midpoint   # "results/sampleid"
outdir    <- opt$outdir

##############
#### Main ####
# sample <- '/mnt/work1/users/pughlab/bin/swgs/results/zygosity/AD/aggregate_filt.csv'
# sample <- '/mnt/work1/users/pughlab/bin/swgs/results/sampleid/chrM/merge.vcf'

# Similarity Matrices
if(ids=='NA') ids <- NULL
sim_mats <- getSampleSimilarity(sample_matrix = sample, samples=ids, matchmode=mode)

# Outputs
## PDF of matrix
pdf(file.path(outdir, paste0("similarity_", mode, ".pdf")))
plotSampleSimilarity(sim_mat = x$sim, n_mat = if(mode=='autosome') x$het else x$n, 
                     midpoint = midpoint, mid_diag=TRUE)
dev.off()

## Written table of matrix
write.table(round(x$sim,3), file = file.path(outdir, paste0("similarity_", mode, ".tsv")),
            row.names = TRUE, col.names = TRUE, sep = "\t")

write.table(if(mode=='autosome') x$het else x$n, 
            file = file.path(outdir, paste0("n_", mode, ".tsv")),
            row.names = TRUE, col.names = TRUE, sep = "\t")
