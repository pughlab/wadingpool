#' Generates a similarity matrix
#' @description Takes a path to a directory that contains
#' the [CHR]_filt.snps files and processes each to generate
#' a similarity matrix using a given similarityFun.  If no
#' similarityFun is given, it generates one using the default
#' Jaccard metric.
#' 
#' @param similarityFun  Similarity function with two parameters, 'i' and 'j'. 
#' @param rm_nocov If no filtering process took place prior, setting this to TRUE
#' will remove all SNPs where no sample have any coverage (Default=FALSE)
#' @param sample_matrix character: path to csv or vcf file to run similarity match on
#' @param matchmode character: Either 'autosome' or 'chrM' (Default=autosome)
#' @param rem_ref boolean: Removes reference matches in chrM sample matching (0/0 and 
#' 0/0) (Default=TRUE)
#' @param samples character: Comma-separated list of ordered samples for the matrix
#' 
#' @importFrom utils read.csv
#' @importFrom stats na.omit
#' @importFrom assertthat assert_that
#' @return
#' Returns a list of similarity matrices, matrix of number of SNPs found between
#' two samples, and a matrix of number of heterozygous SNPs between two samples
#' @export
getSampleSimilarity <- function(sample_matrix, samples, matchmode='autosome',
                       similarityFun=NULL, rm_nocov=FALSE, rem_ref=TRUE){
  assert_that(file.exists(sample_matrix), msg="sample_matrix must be a file path containing the samples to match")
  
  ## Read in the VCF or CSV file
  mat <- switch(matchmode,
                "autosome"={
                  mat <- read.csv(sample_matrix, header = TRUE)
                  if(!is.null(samples)) colnames(mat) <- strsplit(samples, ",")[[1]]
                  mat
                },
                "chrM"={
                  mega_vcf = read.csv(sample_matrix, sep = "\t", comment.char = "#")
                  mat = mega_vcf[,10:ncol(mega_vcf)]
                  if(!is.null(samples)) colnames(mat) <- strsplit(samples, ",")[[1]]
                  mat
                },
                stop("matchmode must be either 'autosome' or 'chrM'"))
  
  ## Run the similarity checks
  mats <- switch(matchmode,
                 "autosome"=.autosomeMatch(mat, rm_nocov, similarityFun),
                 "chrM"=.chrmMatch(mat, rem_ref, similarityFun),
                 stop("matchmode must be either 'autosome' or 'chrM'"))
    
  return(mats)
}

.autosomeMatch <- function(mat, rm_nocov, similarityFun){
  mat <- as.matrix(mat)
  storage.mode(mat) <- 'integer'
  mat <- as.data.frame(mat)
  
  if(rm_nocov){
    # Removes SNPs where no samples have coverage for it
    mat[mat==0] <- NA         # Remove no coverage SNPs
    nonhet_idx <- apply(mat,1,function(i) {
      i <- na.omit(as.integer(i))
      any(i==2) | (any(i==1) & any(i==3))
    })
    mat[!nonhet_idx,] <- NA    # Removes SNPs that have no het cov
  }
  
  if(is.null(similarityFun)){
    similarityFun <- function(i,j){
      # Heterozygous SNPs in i and j
      # Divided by
      # SNPs heterozygous in i or j + 
      #   SNPs that are ref in i and alt in j (vice versa)
      ihet <- i==2
      jhet <- j==2
      #ihom <- (i==1 | i==3)
      # i[which(ihom)] != j[which(ihom)]
      
      sum(ihet & jhet, na.rm=T) / 
        (sum(ihet | jhet, na.rm=T)) #+ sum(i[which(ihom)] != j[which(ihom)], na.rm=T))
    }
  }
  
  getN <- function(i,j){ sum(!is.na(i) & !is.na(j)) }
  getHet <- function(i,j){ sum((i==2) & (j==2),na.rm=T) }
  
  simmat <- allbyall(mat, margin=2, fun=similarityFun)
  nmat <- allbyall(mat, margin=2, fun=getN)
  hetmat <- allbyall(mat, margin=2, fun=getHet)
  
  list("sim"=simmat,
       "n"=nmat,
       "het"=hetmat)
}

.chrmMatch <- function(mat, rem_ref, similarityFun){
  # Parse haplotype caller into the first column (i.e. 0/0, or 1/1, or 0/1).
  sample_geno <- apply(mat, 2, .getGT)
  sample_count <- .cleanGenotype(as.matrix(sample_geno))
  colnames(sample_count) <- colnames(mat)
  
  ## Read all the samples and count the genotypes####
  if(is.null(similarityFun)){
    similarityFun <- function(i,j){
      # Checks that each samples SNP matches the other samples SNP
      comp_idx <- c(1:length(i))
      if(rem_ref){
        #print("Removing SNPs which are ref in both samples")
        refs <- i == '0/0' & j=='0/0'
        comp_idx <- comp_idx[-which(refs)]
      }
      m <- sum(i[comp_idx] == j[comp_idx])
      n <- length(comp_idx)
      jacc <- sum(m)/n
      jacc
      # return(c(jacc,n))
    }
  }
  getN <- function(i,j){ 
    comp_idx <- c(1:length(i))
    if(rem_ref){
      #print("Removing SNPs which are ref in both samples")
      refs <- i == '0/0' & j=='0/0'
      comp_idx <- comp_idx[-which(refs)]
    } 
    return(length(comp_idx))
  }

  jacc_mat  <- allbyall(sample_count, margin=2, fun=similarityFun)
  n_mat     <- allbyall(sample_count, margin=2, fun=getN)
  
  list("sim"=jacc_mat,
       "n"=n_mat)
}

.getGT <- function(x, gt_idx=1){
  split_gt <- strsplit(as.vector(x), split = ":")
  sapply(split_gt, function(i) i[[gt_idx]])
}

.cleanGenotype <- function(b) {
  # GT genotype, encoded as alleles values separated by either of ”/” or “|”,
  # e.g. The allele values are 0 for the reference allele (what is in the
  # reference sequence), 1 for the first allele listed in ALT, 2 for the
  # second allele list in ALT and so on. For diploid calls examples could be
  # 0/1 or 1|0 etc. For haploid calls, e.g. on Y, male X, mitochondrion,
  # only one allele value should be given. All samples must have GT call
  # information; if a call cannot be made for a sample at a given locus, ”.”
  # must be specified for each missing allele in the GT field (for example
  # ./. for a diploid). The meanings of the separators are:
  stopifnot(any(class(b)=='matrix'))
  
  b <- gsub("\\|", "/", b)                 # convert phased to unphased
  b[b == '.'] <- "0/0"                     # convert No-call to HOMREF
  
  # Deeper sanity checks
  splgt <- sapply(strsplit(as.vector(b), "/"), function(i) i[1:2])
  class(splgt) <- 'integer'
  
  # CHECK_1: Flips 1/0 to 0/1
  hetidx <- which(splgt[2,] != splgt[1,])
  rev_hetidx <- splgt[1,hetidx] > splgt[2,hetidx]
  if(any(rev_hetidx)){
    for(i in which(rev_hetidx)){
      splgt[,hetidx[i]] <- rev(splgt[,hetidx[i]])
    }
  }
  
  # CHECK_2: ...
  if(FALSE){
    
  }
  
  genotype <- matrix(apply(splgt, 2, paste, collapse="/"), ncol=ncol(b))
  return(genotype)
}