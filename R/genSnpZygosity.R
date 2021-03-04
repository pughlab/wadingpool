#' Title
#' @description 
#' 
#' @return
#' @export
#'
#' @examples
getSnpZygosity <- function(sampleid, chrs = paste0("chr", c(1:22, "X", "Y"))){
  sampleid <- 'NET-2-001b_02_T_DNA.processed'
  
  sapply(chrs, function(chr){
    vcf_file <- paste0(chr, sampleid, ".vcf")
    
  })
  
  
  self_pkg_name <- 'WadingPool'
  bin_dir <- system.file("bin/", package = self_pkg_name)
  conv.cmd <- paste0(bin_dir, "/categorizeAD.sh", 
                     " -v ", vcf_file, 
                     " -o ", out_file)
  conv.cmd.res <- try(system(command = conv.cmd, intern = TRUE))
}
