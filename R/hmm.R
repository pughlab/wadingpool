#' Fit an HMM
#' @description Runs depmixs4 to fit an HMM to a set of observations
#' given an expected number of states. It will then configure the 
#' emission dataframe to be compatible with downstream analysis
#'
#' @param sample String containing sample ID in columns
#' @param het_cnt Dataframe for sample observations
#' @param states Integer for expected number of states
#' @param family Family of distribution to use (default=poisson())
#' @param tiles GenomicRanges object for corresponding chrs and bins.  If none is given, it
#' will fill in the genomic poisition randomly (Default=NULL)
#' @param is.sim Boolean flag for whether data is simulation or not (Default=FALSE)
#'
#' @importFrom depmixS4 depmix
#' @importFrom depmixS4 fit
#' @importFrom depmixS4 posterior
#' @importFrom GenomeInfoDb seqnames
#' 
#' @return
#' Dataframe containing the emission proabbility for each state and a fitted HMM
#' @export
fitHMM <- function(sample, het_cnt, states=2, family=poisson(), tiles=NULL, is.sim=FALSE){
  ## Fit HMM with best estimate of k
  sample <- 'obs'
  mod <- depmix(as.formula(paste0(sample, '~ 1')), 
                data = as.data.frame(het_cnt), nstates=states, family=family)
  fit.mod <- fit(mod)
  
  # Label the estimated States based on posterior prob
  est.states     <- cbind("obs"=het_cnt[,sample], posterior(fit.mod))
  state_obs      <- split(est.states$obs, est.states$state)
  state_obs_mean <- sort(sapply(state_obs, mean))
  state_obs_mean <- setNames(c("LOH", paste0("H", c(1:(length(state_obs_mean)-1)))), names(state_obs_mean))
  state_hash     <- as.list(state_obs_mean)
  
  est.states$state.label <- unlist(state_hash[as.character(est.states$state)])
  # est.states$state.label <- rep('Het', nrow(est.states))
  # est.states$state.label[est.states$state == which.min(state_obs_mean)] <- 'LOH'
  
  # Label the chromosomes
  if(!is.null(tiles)){
    est.states$chrom <- as.character(seqnames(tiles)) #as.character(seqnames(tiles[as.integer(names(ov_spl)),]))
    est.states$bin <- names(tiles) #names(tiles[as.integer(names(ov_spl)),])
  } else {
    est.states$chrom <- paste0("chrom", sort(rep(1:22, 200)))[1:nrow(est.states)]
    est.states$bin <- paste0("bin", c(1:nrow(est.states)))
  }
  
  if(is.sim){
    est.states$act.state.label <- het_cnt$state
    est.states$act.state <- 0
    est.states <- .genStateMap(est.states, est_cols = c('state', 'state.label'), 
                               act_cols = c('act.state', 'act.state.label'))
  }
  
  return(est.states$model)
}
