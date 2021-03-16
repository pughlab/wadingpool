devtools::install_github("quevedor2/WadingPool", ref = "dev")

library(WadingPool)
PDIR <- '/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/gatk_common/output'
setwd(PDIR)

## Run Heterozygous Count Simulation
grch38_dbsnp <- file.path('/mnt/work1/users/home2/quever/pughlab/references/dbsnp/dbsnp_b151_GRCh38p7/chr_bed',
                          'all.common_all_20180418.bed')
exp_hets <- calcExpectedHets(grch38_dbsnp, coverage=seq(0,3,by=0.1))
samp1 <- getExpectedN(mu = exp_hets['mean',], se = exp_hets['se',], n = exp_hets['n',])
samp2 <- getExpectedN(mu = exp_hets['meansq',], se = exp_hets['se',], n = exp_hets['n',])
boxplot(t(samp2[-2,4:6])/1e3, log = 'y', las=2)


library(nlstools)
nlsdata <- data.frame("cov"=as.numeric(rownames(samp2)),
                      "SNPs"=samp2$n_mean)[-1,]
nlsfitSS <- nls(SNPs ~ SSasymp(cov, Asym, R0, lrc),
                data=nlsdata)

pred_coverage = seq(0,10,by=0.01)
pred_snp <- predict(nlsfitSS,list(cov=pred_coverage))
Asym_coef<-summary(nlsfitSS)$coefficients[1] #100% saturation
R0_coef<-summary(nlsfitSS)$coefficients[2]
lrc_coef<-summary(nlsfitSS)$coefficients[3]
sat95<-Asym_coef*0.95 #peaks at 95% saturation
sat99<-Asym_coef*0.99 #peaks at 99% saturation
SN_sat95<--(log((sat95-Asym_coef)/(R0_coef-Asym_coef))/exp(lrc_coef)) # Coverage to reach 95% saturation
SN_sat99<--(log((sat99-Asym_coef)/(R0_coef-Asym_coef))/exp(lrc_coef)) # Coverage to reach 99% saturation






suffix <- 'allelicCounts.tsv'
files <- list.files(PDIR, pattern=paste0(".", suffix, "$"))
chrs <- 'all'
samples <- gsub(paste0("^", unique(gsub("[0-9XY]", "", chrs)), 
                       ".*?\\.(.*)",
                       paste0(".", suffix, "$")), "\\1", files)
# samples <- unique(gsub("^chr.*?\\.(.*).vcf$", "\\1", files))
# chrs <- paste0("chr", c(1:22, "X", "Y"))

snp_results <- lapply(samples, genSnpZygosity, caller='GATK', chrs=chrs)

s_paths <- file.path("tmp", gsub("[^a-zA-Z0-9]", "", samples), 
                     paste0("[ID].", samples, "_out.tsv"))
agg_results <- aggregateAndFilter(s_paths, num_samples=2, chrs=chrs,
                                  dbsnp_file='all.common_all_20151104.bed',
                                  dbsnp_dir='/mnt/work1/users/pughlab/references/dbsnp/dbsnp_b146_GRCh37p13/chr_bed') 

filt_dir <- file.path("tmp", "combined", "filt")
getSampleSimilarity(filt_dir, samples)
