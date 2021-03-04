#' Simplifies the SNP Zygosity from MuTect
#' 
#' @description Takes a single sample Mutect VCF that has force called all dbSNP
#' sites and reduces it to a single column of 1,2,3 values corresponding to 
#' REF.HOMOZYGOUS, HETEROZYGOUS, and ALT.HOMOZYGOUS respectively using Perl. 
#' If multiple samples are given, it will concatenate the samples into a single matrix.
#' 
#' @param sampleid Vector of sample IDs for the filename [chr].[sampleid].vcf
#' @param chrs Chrs for the filename [chr].[sampleid].vcf
#' 
#' @return
#' Matrix of parsed output values
#' @export
genSnpZygosity <- function(sampleid, chrs = paste0("chr", c(1:22, "X", "Y"))){
  # sampleid <- 'NET-2-001b_02_T_DNA.processed'
  sampleid_dir <- file.path("tmp", gsub("[^a-zA-Z0-9]", "", sampleid))
  dir.create(sampleid_dir, recursive = TRUE, showWarnings = FALSE)
  
  self_pkg_name <- 'WadingPool'
  bin_dir <- system.file("bin/", package = self_pkg_name)
  
  categorize_results <- sapply(chrs, function(chr){
    vcf_file <- paste0(chr, ".", sampleid, ".vcf")
    out_file <- file.path(sampleid_dir, paste0(chr, ".", sampleid, "_out.tsv"))
    
    ## Executes the categorizeAD.sh shell script that parses the single-sample
    #  MuTect VCF output 
    categorize_cmd <- paste0("sh ", bin_dir, "/categorizeAD.sh", 
                             " -v ", vcf_file, 
                             " -o ", out_file)
    categorize_res <- try(system(command = categorize_cmd, intern = TRUE))
    
    return(c("Time"=as.numeric(parseOutput(categorize_res, "Runtime")),
             "In"=as.character(parseOutput(categorize_res, "vcf")),
             "Out"=as.character(parseOutput(categorize_res, "out"))))
  })
  categorize_results <- t(categorize_results)
  
}
