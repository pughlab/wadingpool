.check <- function(x, check, ...){
  switch(check,
         plotHMM.est_cols=.chk_est_cols(x, ...))
}

#########################
#### hmm_viz.PlotHMM ####
.chk_est_cols <- function(x, m){
  assert_that(length(x)==2, is.character(x), all(x %in% colnames(m)),
              msg="est_cols must be a 2 element character vector for col names in model")
  assert_that(is.numeric(m[,x[1]]), 
              msg="est_cols[1] must correspond to a numeric state c(1,2,3,...,n)")
  assert_that(is.character(m[,x[2]]), 
              msg="est_cols[2] must correspond to a character label for each state")
}

#################################
#### dataProcessing.parseSeg ####
.chk_seg <- function(seg, st_end, event){
  
}