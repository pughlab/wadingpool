#' Aggregate samples and filter SNPs
#' 
#' @description Works in conjunction with the genSnpZygosity() function.
#' This function will take the output from that function and concatenate
#' each samples SNP results into a csv. It will then go through 
#' each SNP to identify whether X or more samples have coverage on that
#' SNP and isolate the genomic location and genotypes for those SNPs
#' 
#' @param s_paths Sample SNP paths from genSnpZygosity() 
#' @param dbsnp_dir Directory containing the chr-dbSNP bed files
#' @param num_samples Number of samples required to have that SNP [Default=2]
#' @param dbsnp_file Main dbSNP file before splitting [default=all.common_all_20151104.bed]
#' @param chrs List of chrs [Default=chr1, chr2, .., chrY]
#'
#' @return
#' A vector of the input, output files and runtime
#' @export
aggregateAndFilter <- function(s_paths, dbsnp_dir, num_samples=2, 
                               dbsnp_file='all.common_all_20151104.bed',
                               chrs = paste0("chr", c(1:22, "X", "Y"))){
  
  # dbsnp_dir <- '/mnt/work1/users/pughlab/references/dbsnp/dbsnp_b146_GRCh37p13/common/bed'
  # dbsnp_file <- 'all.common_all_20151104.bed'
  # num_samples <- 1
  
  self_pkg_name <- 'WadingPool'
  bin_dir <- system.file("bin/", package = self_pkg_name)
  
  # Set Paths
  comb_path <- file.path("tmp", "combined")
  dir.create(file.path(comb_path, "idx"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(comb_path, "filt"), recursive = TRUE, showWarnings = FALSE)
  
  aggregate_res <- sapply(chrs, function(chr){
    # Set I/O filenames
    agg_out <- paste0(chr, ".snps")         # combined all SNPs
    line_snps <- paste0(chr, "_lines.txt")  # Line #s for filtered SNPs
    filt_snps <- paste0(chr, "_filt.snps")  # Filtered SNPs
    pos_snps <- paste0(chr, "_pos.snps")    # Genomic loci of SNPs
    dbsnpf <- gsub("^all", chr, dbsnp_file) # Chr specific dbSNP reference file
    
    start_time <- Sys.time()                # Start the timer
    
    ## Pastes individual SNP columns together into a csv
    # paste -d"," [chrX.file1.tsv] [chrX.file1.tsv] > chrX.snps
    chr_s_paths <- gsub("\\[ID\\]", chr, s_paths)
    paste_cmd <- paste0('paste -d "," ',
                        paste(chr_s_paths, collapse=" "),
                        " > ", file.path(comb_path, agg_out))
    paste_res <- try(system(command = paste_cmd, intern = TRUE))
    
    ## Identify SNPs that are covered in X number of samples
    awk_cmd <- paste0('awk -F\',\' \'{ for(i=1; i<=NF;i++) if ($i!=0) j++; ', 
                      'if (j >= ', num_samples, ') print NR; j=0 }\' ',
                      file.path(comb_path, agg_out), " > ", 
                      file.path(comb_path, "idx", line_snps))
    awk_res <- try(system(command = awk_cmd, intern = TRUE))
    
    ## Isolate the SNPs that are covered in X number of samples
    # perl getLinesFromFile.py ${linesnps} ${snps} > ${filtsnps}
    extlines_cmd <- paste0("perl ", bin_dir, "/getLinesFromFile.py ", 
                           file.path(comb_path, "idx", line_snps), " ", 
                           file.path(comb_path, agg_out), 
                           "> ", file.path(comb_path, "filt", filt_snps))
    extlines_res <- try(system(command = extlines_cmd, intern = TRUE))
    
    ## Isolate the genomic loci for the SNPs from dbSNP reference:
    # perl getLinesFromFile.py ${linesnps} ${dbsnpdir}/${dbsnpfile} > ${possnps}
    dbsnpiso_cmd <- paste0("perl ", bin_dir, "/getLinesFromFile.py ", 
                           file.path(comb_path, "idx", line_snps), " ", 
                           file.path(dbsnp_dir, dbsnpf), 
                           "> ", file.path(comb_path, "filt", pos_snps))
    dbsnpiso_res <- try(system(command = dbsnpiso_cmd, intern = TRUE))
    
    end_time <- Sys.time()                  # End the timer
    runtime <- as.integer(difftime(end_time, start_time, units='secs'))
    
    
    return(c("Time"=runtime,
             "In"=paste(chr_s_paths, collapse=","),
             "Out"=paste(c(file.path(comb_path, agg_out), 
                           file.path(comb_path, "filt", filt_snps), 
                           file.path(comb_path, "filt", pos_snps)), 
                         collapse=",")))
  })
  return(aggregate_res)
}

