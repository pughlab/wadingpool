# devtools::install_github("quevedor2/WadingPool", ref = "dev")
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(WadingPool))
suppressPackageStartupMessages(library(GenomicRanges))

####################
#### Parameters ####
option_list <- list(
  make_option(c("-p", "--hetposf"), type="character", default='results/zygosity/AD/aggregate_pos.txt',
              help="File containing genomic position of heterozygous SNPs [default=%default]"),
  make_option(c("-f", "--hetf"), type="character", default='results/zygosity/AD/aggregate_filt.csv',
              help="File containing counts for heterozygous SNPs [default=%default]"),
  make_option(c("-c", "--cnpath"), type="character", default='results/cnv/ichorcna',
              help="Path to directory containing copy-number seg files [default=%default]"),
  make_option(c("-o", "--outdir"), type="character", default='results/zygosity/wadingpool',
              help="Path to out directory [default=%default]"),
  make_option(c("-g", "--genome"), type="character", default='hg19',
              help="Genome build to use: hg19, hg38, mm10 [default=%default]"),
  make_option(c("-m", "--maxstate"), type="integer", default=300,
              help="Maximum number of segments for zygosity HMM before reducing number of states [default=%default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ucsc_dir <- '/mnt/work1/users/pughlab/references/ucsc/hg19'
# pdir     <- '/mnt/work1/users/pughlab/bin/swgs/'
hetcnt_pos_f  <- opt$hetposf  # file.path(pdir, "results/zygosity/AD", "aggregate_pos.txt")
hetcnt_f      <- opt$hetf     # file.path(pdir, "results/zygosity/AD", "aggregate_filt.csv")
cn_path       <- opt$cnpath   # "results/cnv/ichorcna"
outdir        <- opt$outdir

model    <- opt$genome        # 'hg19'
maxstate <- opt$maxstate      # 300
chrs     <- paste0("chr", c(1:22,"X", "Y"))

########################
#### Bin the genome ####
if(model=='hg19'){
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome <- BSgenome.Hsapiens.UCSC.hg19
} else if(model=='hg38'){
  stop("Not tested")
  library(BSgenome.Hsapiens.UCSC.hg38)
  genome <- BSgenome.Hsapiens.UCSC.hg38
} else if(model=='mm10'){
  stop("Not tested")
  library(BSgenome.Mmusculus.UCSC.mm10)
  genome <- BSgenome.Mmusculus.UCSC.mm10
}
tiles <- tileGenome(seqinfo(genome)[chrs], tilewidth=1e6, 
                    cut.last.tile.in.chrom=TRUE)
names(tiles) <- paste0("bin", c(1:length(tiles)))


######################################
#### Read in Het-Cnt and Bin data ####
# Read in het-cnt positional data
pos <- read.table(hetcnt_pos_f, sep="\t", header=FALSE)
colnames(pos) <- c('chrom', 'start', 'end')
grpos <- makeGRangesFromDataFrame(pos)
rm(pos)

# Read in Het-Cnt data
ad <- read.csv(hetcnt_f, header=TRUE, check.names = FALSE)

# Bin the Het-Cnt data
ov <- findOverlaps(tiles, grpos)
ov_spl <- split(ov, queryHits(ov))
tiles_ov <- tiles[unique(sort(queryHits(ov)))]
het_cnt <- t(sapply(ov_spl, function(i) colSums(ad[subjectHits(i),]==2)))
colnames(het_cnt) <- colnames(ad)

# Sanitize data
het_cnt[is.na(het_cnt)] <- 0
het_cnt <- apply(het_cnt,2,sanitizeData)


##############################
#### Read in CN Seg files ####
cndat <- lapply(setNames(colnames(het_cnt), colnames(het_cnt)), function(sample){
  # cn_path <- "results/cnv/ichorcna"
  seg_id <- paste0(sample, ".seg")
  if(file.exists(file.path(cn_path, seg_id))){
    segf <- file.path(cn_path, seg_id)
  } else if(file.exists(file.path(cn_path, sample, seg_id))) {
    # ichorCNA path
    segf <- file.path(cn_path, sample, seg_id)
  } else {
    stop(paste0(seg_id, " could not be found"))
  }

  seg_state <- parseSeg(segf, get_seg = TRUE)
  seg_gr <- makeGRangesFromDataFrame(seg_state$seg, keep.extra.columns = T)
  seqlevelsStyle(seg_gr) <- 'UCSC'
  list("state"=seg_state$states, "seg"=seg_gr)
})
states <- sapply(cndat, function(i) i$state)
segs <- as(lapply(cndat, function(i) i$seg), "GRangesList")

#################################
#### Fit HMM and plot/report ####
hmm_fit <- lapply(setNames(colnames(het_cnt), colnames(het_cnt)), function(sample){
  state <- states[[sample]] 
  
  # Model fitting (iterate until number of states < max number of states)
  numstates <- maxstate
  while(numstates >= maxstate){
    hmm1 <- WadingPool::fitHMM(sample=sample, het_cnt = as.data.frame(het_cnt), states=state, family=gaussian(), 
                               tiles=tiles_ov, is.sim = FALSE, ret.raw=TRUE, multi=FALSE, 
                               tr_same = c(0.9, 0.99), tr_switch = c(0.001, 0.1))
    print(paste0(sample, ": (", state, ") - ", length(rle(hmm1$df$state)$lengths)))
    state <- states[[sample]]  - 1
    numstates <- length(rle(hmm1$df$state)$lengths)
  }
  
  #OUTPUTS
  ## CN and Zyg-HMM visualizations
  pdf(file.path(outdir, paste0("hmmfit_", sample, ".pdf")))
  seg <- plotCnZyg(segs[[sample]], combZygPos(zyg=hmm1$df, zygpos=tiles_ov), tiles, main=sample)
  WadingPool::plotHMM(model=hmm1$df, est_cols = c('state', 'state.label'), xval = 'bin')
  dev.off()
  
  ## Zyg-HMM PP table output
  write.table(hmm1$df, file=file.path(outdir, paste0("hmmfit_", sample, ".tsv")), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  ## Zyg-CN seg output
  write.table(seg, file=file.path(outdir, paste0("cn-zyg_", sample, ".tsv")), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  ## model output
  model <- hmm1$model
  save(model, file=file.path(outdir, paste0("hmmfit_", sample, ".rda")))
  
  hmm1[['seg']] <- seg
  return(hmm1)
})

