devtools::install_github("quevedor2/WadingPool", ref = "dev")

library(WadingPool)
PDIR <- '/mnt/work1/users/pughlab/projects/NET-SEQ/shallow_wgs/variant_calling/gatk_common/output'
setwd(PDIR)

# Simulation
{
## Run Heterozygous Count Simulation
grch38_dbsnp <- file.path('/mnt/work1/users/home2/quever/pughlab/references/dbsnp/dbsnp_b151_GRCh38p7/chr_bed',
                          'all.common_all_20180418.bed')
load("~/exp_hets.rda")
exp_hets <- calcExpectedHets(grch38_dbsnp, coverage=seq(0,3,by=0.1))
samp1 <- getExpectedN(mu = exp_hets['mean',], 
                      se = exp_hets['se',], 
                      n = exp_hets['n',])
samp2 <- getExpectedN(mu = exp_hets['meansq',], 
                      se = exp_hets['se',], 
                      n = exp_hets['n',])


singlesample=data.frame("cov"=as.numeric(rownames(samp2)),
                        "SNPs"=samp1$n_mean)[-1,]
single_sat <- saturationCurve(singlesample,
                              pred = seq(0, 10, by=0.01),
                              S = sort(c(seq(0,1,by=0.1), 0.95, 0.99)),
                              Xin = seq(0, 2, by=0.1))


plot(x = rownames(samp1)[-1], y=samp1[-1,]$n_mean,
     xlim=c(0,5), ylim=c(0,single_sat$Asymp_coef), las=2, 
     ylab='', xlab='Coverage')
axis(side = 4, at = single_sat$Yfit, labels = names(single_sat$Yfit), las=2)
rect(xleft = as.numeric(rownames(samp1)[-1]) - 0.02, ybottom = samp1[-1,]$n_low, 
     xright = as.numeric(rownames(samp1)[-1]) - 0.02, ytop = samp1[-1,]$n_high)
lines(x=as.numeric(names(single_sat$pred)), y=single_sat$pred, lty=2)
lines(x = c(single_sat$Xfit['0.95'], single_sat$Xfit['0.95']), 
      y = c(-100, single_sat$Yfit['0.95']), col='red', lty=2)
lines(x = c(single_sat$Xfit['0.95'], 100), 
      y = c(single_sat$Yfit['0.95'], single_sat$Yfit['0.95']), col='red', lty=2)



singlesample=data.frame("cov"=as.numeric(rownames(samp2)),
                        "SNPs"=samp1$n_mean)[-1,]
single_sat <- saturationCurve(singlesample,
                              pred = 0.4, # Predicting for 1x coverage
                              S = 0.95, # Predicting Num. of HET SNPs at 0.95 saturation
                              Xin = 0.4) # Predicting saturation for 1x coverage

twosamples=data.frame("cov"=as.numeric(rownames(samp2)),
                      "SNPs"=samp2$n_mean)[-1,]
two_sat <- saturationCurve(twosamples,
                           pred = 0.4, # Predicting for 1x coverage
                           S = 0.95, # Predicting Num. of HET SNPs at 0.95 saturation
                           Xin = 0.4) # Predicting saturation for 1x coverage

twosamples=data.frame("cov"=as.numeric(rownames(samp2)),
                      "SNPs"=samp2$n_mean)[-1,]
# exp_cov <- c(seq(0, 0.01, by=0.001), seq(0.01, 0.1, by=0.01), seq(0.1, 1, by=0.1))
exp_cov <- seq(0, 2, by=0.1)
two_sat <- saturationCurve(twosamples,
                           pred = seq(0, 10, by=0.01),
                           S = sort(c(seq(0,1,by=0.1), 0.95, 0.99)),
                           Xin = exp_cov)

samp <- samp2
sat <- two_sat
par(mar=c(5.1, 6, 4.1, 4.1))
plot(x = rownames(samp)[-1], y=log10(samp[-1,]$n_mean),
     xlim=c(0,5), ylim=c(0,log10(sat$Asymp_coef)), las=2, 
     ylab='', xaxt='n', yaxt='n', xlab='Coverage')
axis(side = 4, at = log10(sat$Yfit), labels = names(sat$Yfit), las=2)
mtext("Saturation", side=4, line=2, cex.lab=1)
axis(side = 1, at = c(seq(0, 1, by=0.1), seq(1, 5, by=1)), 
     labels = c(seq(0, 1, by=0.1), seq(1, 5, by=1)), las=2)
axis(side = 2, at = seq(1:log10(max(samp[-1,]$n_mean))), 
     labels = 10^(seq(1:log10(max(samp[-1,]$n_mean)))), las=2)
mtext("Number of Het. SNPs", side=2, line=3, cex.lab=1)
# rect(xleft = as.numeric(rownames(samp)[-1]) - 0.02, ybottom = samp[-1,]$n_low, 
#      xright = as.numeric(rownames(samp)[-1]) - 0.02, ytop = samp[-1,]$n_high)
lines(x=as.numeric(names(sat$pred)), y=log10(sat$pred), lty=2)
lines(x = c(sat$Xfit['0.95'], sat$Xfit['0.95']), 
      y = c(-100, log10(sat$Yfit['0.95'])), col='red', lty=2)
lines(x = c(sat$Xfit['0.95'], 100), 
      y = c(log10(sat$Yfit['0.95']), log10(sat$Yfit['0.95'])), col='red', lty=2)


two_sat$S[two_sat$S < 0] <- 0
two_sat_response <- saturationCurve(twosamples,
                           pred = two_sat$S)
two_sat_response$pred[two_sat_response$pred<0] <- 0
pred_snps <- data.frame("Coverage"=exp_cov, 
                        "Saturation"=two_sat$S,
                        "Num.of.SNPs"=two_sat_response$pred)
}



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

sample <- '/mnt/work1/users/pughlab/bin/swgs/results/zygosity/AD/aggregate_filt.csv'
x <- getSampleSimilarity(sample_matrix = sample, matchmode='autosome')

sample <- '/mnt/work1/users/pughlab/bin/swgs/results/sampleid/chrM/merge.vcf'
x <- getSampleSimilarity(sample_matrix = sample, 
                         samples=c('net-001a,net-001b,net-002,net-037'), 
                         matchmode='chrM')