#' Calculates expected number of heterozygous SNPs
#' @description Takes a dbSNP file that is converted into the .bed format
#' by the setup_dbSNP.sh script and runs simulations to calculate the 
#' expected number of heterozygous SNPs in a single sample given a target
#' coverage fitted to a Poisson distribution.
#'
#' @param dbsnp_file Path to dbSNP bed file (from setup_dbSNP.sh)
#' @param coverage Target coverage(s) and lambda param of rpois()
#' @param max_cnt Count to set all dbSNP at for max point (default=2)
#' @importFrom stats rpois
#' @importFrom stats setNames
#' @importFrom utils read.table
#' @return
#' A matrix for the input coverages and the mean, mean^2 (for 2 samples 
#' overlapping), standard error heterozygous SNPs (Expected values) and
#' the number of SNPs
#' @export
calcExpectedHets <- function(dbsnp_file, coverage=c(1), max_cnt=2){
  # setwd('/mnt/work1/users/home2/quever/pughlab/references/dbsnp/dbsnp_b151_GRCh38p7/chr_bed')
  # dbsnp_file='all.common_all_20180418.bed'
  
  dbsnp <- read.table(dbsnp_file, header = F,
                      stringsAsFactors = F, check.names = F)
  colnames(dbsnp)[c(1:5)] <- c('chrom', 'start', 'end', 'rs', 'caf')
  dbsnp$adjcaf <- (0.5 - abs(0.5 - dbsnp$caf))    # Probabilty of Het
  
  ## Function to calculate the probabilty of getting a heterozygous
  # SNP for a given coverage using a geometric series; for example:
  #   - AF of 0.5 with Cov of 2, P(Het)=0.5
  #   - AF of 0.5 with Cov of 3, P(Het)=0.75
  #   - AF of 0.5 with Cov of 4, P(Het)=0.875
  probHet <- function(p, n){
    x <- matrix(c(p, n), ncol=2, byrow = FALSE)
    phet <- apply(x, 1, function(i) round(sum(i[1]^(1:(i[2]-1))),3))
    phet[phet>1] <- 0
    return(phet)
  }
  
  coverage <- c(-1,coverage)
  exp_hets <- sapply(setNames(coverage,coverage), function(cov){
    if(cov==-1){
      print(paste0("Setting max where each SNP is covered by ", max_cnt, " reads"))
      dbsnp$cov <- max_cnt
    } else {
      print(paste0("Coverage: ", cov))
      # Simulate coverage using a poisson distribution
      dbsnp$cov <- rpois(n = nrow(dbsnp), lambda = cov)
    }
    
    # Use a geometric series to calculate the P(Het) per SNP
    dbsnp$p_het <- probHet(dbsnp$caf, dbsnp$cov)
    
    ## Calculate expected value of heterozygous SNPs and 
    # the 95% confidence itnerval around it
    exphet_mean <- mean(dbsnp$p_het)
    exphet_se <- sqrt(exphet_mean/nrow(dbsnp))
    
    c("mean"=exphet_mean,
      "meansq"=exphet_mean^2,
      "se"=exphet_se,
      'n'=nrow(dbsnp))
  })
  return(exp_hets)
}

#' Calculates the expected number of SNPs 
#' @description Calculates the 95% confidence interval for the number of SNPs
#' given a total N, as well as the distance between these SNPs
#' 
#' @param mu Mean
#' @param se Standard error
#' @param n N
#' @param genome_size Default (3.2*10^9)
#' @importFrom stats setNames
#' @return 
#' 95% confidence for number of SNPs given mu, se, and number of SNPs
#' 95% confidence for genomic distance between those SNPs
#' @export
getExpectedN <- function(mu, se, n, genome_size=3.2*10^9){
  nval <- c("n_low"=(mu - (1.96*se))*n,
            "n_mean"=mu*n,
            "n_high"=(mu + (1.96*se))*n)
  return(c(nval, 
           setNames((genome_size/nval), c("dist_low", "dist_mean", "dist_high"))))
}
