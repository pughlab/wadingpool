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
#' @param ret.raw Boolean flag to return the fitted model depmixS4 object
#'
#' @importFrom depmixS4 depmix
#' @importFrom depmixS4 fit
#' @importFrom depmixS4 posterior
#' @importFrom depmixS4 em.control
#' @importFrom GenomeInfoDb seqnames
#' 
#' @return
#' Dataframe containing the emission proabbility for each state and a fitted HMM
#' @export
fitHMM <- function(sample, het_cnt, states=2, family=poisson(), tiles=NULL, is.sim=FALSE, ret.raw=FALSE){
  ## Fit HMM with best estimate of k
  # sample <- 'obs'
  mod <- depmix(as.formula(paste0(sample, '~ 1')), data = as.data.frame(het_cnt), 
                nstates=states, family=family, trstart = runif(states^2))
  fit.mod <- fit(mod, emcon=em.control(random.start=F))
  
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
    est.states$act.state <- NA
    est.states <- .genStateMap(est.states, est_cols = c('state', 'state.label'), 
                               act_cols = c('act.state', 'act.state.label'))
  }
  
  ret_obj <- est.states$model
  if(ret.raw){
    ret_obj <- fit.mod
  }
  return(ret_obj)
}


#' @description supporting function to get the mapping between
#' state and state labels
#'
#' @param model from fitHmm()
#' @param est_cols 2 element character vector for the estimated [state, state.label] columns
#' @param act_cols If known, 2 element character vector for the actual [state, state.label] columns
#'
#' @return
#' 2 element list:
#'   'model': Model with the actual states mapped to be concordant with estimated labels
#'   'map': Unique mapping between state and state-labels
.genStateMap <- function(model, est_cols, act_cols=NULL){
  uniq_map <- unique(model[,est_cols])
  uniq_map <- uniq_map[order(uniq_map[,est_cols[2]]),]
  
  if(!is.null(act_cols)){
    # Re-assign the "Actual" state - state_label mapping to be concordant with the 
    # "Estimated" state - state_label mapping
    .mapStates <- function(x, uniq_map){
      uniq_map[match(x, uniq_map[,2]),1]
    }
    
    ## Map state value for estimate state-labels to actual state-labels
    model[,act_cols[1]] <- model[,act_cols[2]] %>%
      map(.mapStates, uniq_map) %>% 
      unlist
    
    ## Check if all actual states have a proper mapping
    uniq_map_act <- unique(model[,act_cols])
    uniq_map_act <- uniq_map_act[order(uniq_map_act[,act_cols[2]]),]
    na_idx <- is.na(uniq_map_act[,1])
    if(any(na_idx)){
      # Case if actual states are more than the estimated states
      uniq_map_act[which(na_idx),1] <- c(1:10)[-sort(uniq_map_act[,1])][1:sum(na_idx)]
      model[,act_cols[1]] <- model[,act_cols[2]] %>%
        map(.mapStates, uniq_map_act) %>% 
        unlist
    }

    if(!all(unique(model[,act_cols[2]]) %in% uniq_map[,2])) warning("Actual labels do not match the estimated labels")
    assert_that(!any(is.na(uniq_map_act[,1])), !any(is.na(uniq_map[,1])),
                msg="Unmapped states in estimated or actual states")
  }
  return(list("map"=uniq_map, "model"=model, "act.map"=uniq_map_act))
}
