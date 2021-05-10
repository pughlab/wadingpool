#' Generates a similarity matrix
#' @description Takes a path to a directory that contains
#' the [CHR]_filt.snps files and processes each to generate
#' a similarity matrix using a given similarityFun.  If no
#' similarityFun is given, it generates one using the default
#' Jaccard metric.
#' 
#' @param filt_dir  Path to directory containing [CHR]_filt.snps
#' @param samples  Sample vector corresponding to the column names of _filt.snps
#' @param similarityFun  Similarity function with two parameters, 'i' and 'j'. 
#' @param pattern Pattern for the aggregate SNPs file (Default="filt.snps$")
#' @param rm_nocov If no filtering process took place prior, setting this to TRUE
#' will remove all SNPs where no sample have any coverage (Default=FALSE)
#' If no similarity function is given, it default to a Jaccard metric
#' @importFrom utils read.csv
#' @importFrom stats na.omit
#' @return
#' Returns a list of similarity matrices, matrix of number of SNPs found between
#' two samples, and a matrix of number of heterozygous SNPs between two samples
#' @export
getSampleSimilarity <- function(filt_dir, samples,
                       similarityFun=NULL, pattern="filt.snps$", rm_nocov=FALSE){
  # similarityFun=NULL
  if(is.null(filt_dir)) filt_dir <- getwd()
  filt_files <- list.files(filt_dir, pattern=pattern)
  
  jaccmats <- lapply(filt_files, function(f){
    # f <- filt_files[1]
    mat <- read.csv(file.path(filt_dir, f), header = FALSE)
    mat <- apply(mat,2,as.integer)
    colnames(mat) <- samples
    
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
  })
  return(jaccmats)
}