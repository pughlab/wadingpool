#' Title
#'
#' @param filt_dir  to fill in
#' @param samples  to fill in
#' @param similarityFun  to fill in
#' @importFrom utils read.csv
#' @importFrom stats na.omit
#' @return
#' to fill in
#' @export
getSampleSimilarity <- function(filt_dir, samples,
                       similarityFun=NULL){
  # similarityFun=NULL
  filt_files <- list.files(filt_dir, pattern="filt.snps$")
  
  jaccmats <- lapply(filt_files, function(f){
    # f <- filt_files[1]
    mat <- read.csv(file.path(filt_dir, f), header = FALSE)
    colnames(mat) <- samples
    mat[mat==0] <- NA         # Remove no coverage SNPs
    nonhet_idx <- apply(mat,1,function(i) {
      i <- na.omit(as.integer(i))
      any(i==2) | (any(i==1) & any(i==3))
    })
    mat[!nonhet_idx,] <- NA    # Removes SNPs that have no het cov
    
    if(is.null(similarityFun)){
      similarityFun <- function(i,j){
        # Heterozygous SNPs in i and j
        # Divided by
        # SNPs heterozygous in i or j + 
        #   SNPs that are ref in i and alt in j (vice versa)
        ihet <- i==2
        jhet <- j==2
        ihom <- (i==1 | i==3)
        i[which(ihom)] != j[which(ihom)]
        
        sum(ihet & jhet, na.rm=T) / 
          (sum(ihet | jhet, na.rm=T)) #+ sum(i[which(ihom)] != j[which(ihom)], na.rm=T))
      }
    }
    
    getN <- function(i,j){ sum(!is.na(i) & !is.na(j)) }
    
    simmat <- allbyall(mat, margin=2, fun=similarityFun)
    nmat <- allbyall(mat, margin=2, fun=getN)
    
    list("jacc"=simmat,
         "n"=nmat)
  })
}