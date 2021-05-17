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
#' @param nstart Number of starts to try out for EM (Default=10)
#' @param multi Boolean whether to use depmixs4::multistart for fitting model
#'
#' @importFrom depmixS4 depmix
#' @importFrom depmixS4 fit
#' @importFrom depmixS4 posterior
#' @importFrom depmixS4 em.control
#' @importFrom depmixS4 multistart
#' @importFrom GenomeInfoDb seqnames
#' 
#' @return
#' Dataframe containing the emission proabbility for each state and a fitted HMM
#' @export
fitHMM <- function(sample, het_cnt, states=2, family=poisson(), 
                   nstart=10, tiles=NULL, is.sim=FALSE, ret.raw=FALSE,
                   multi=FALSE){
  # Set base transition and emission expectations
  trstart <- matrix(runif(states^2, min=0, max=0.4), ncol=states)
  diag(trstart) <- runif(n = states, min=0.7, max=0.9)  # Non-changing states
  if(family$family == 'gaussian'){
    # gaussians centered around 0, separated by 0.5, sd=0.1
    obs_range   <- round(diff(quantile(het_cnt[,sample], c(0.25, 0.75))),1)
    obs_spread  <- round(obs_range/states,2)
    respstart <- cumsum(rep(obs_spread, states))
    # respstart <- seq(0.5, states*0.5, by=0.5)
    respstart <- as.vector(rbind(respstart - median(respstart), rep(0.1, states)))
  }
  

  mod <- depmix(as.formula(paste0('`', sample, '` ~ 1')), data = as.data.frame(het_cnt), 
                nstates=states, family=family, trstart = unlist(trstart), respstart = unlist(respstart))
  # 
  fit.mod <- tryCatch({
    if(multi){
      multistart(mod, nstart=nstart, initIters=10)
    } else {
      fit(mod, emcon=em.control(random.start=F))
    }
  }, error=function(e){NULL})
  if(is.null(fit.mod)) return(fit.mod)
  
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
    est.states$chrom <- paste0("chr", sort(rep(1:22, 200)))[1:nrow(est.states)]
    est.states$bin <- paste0("bin", c(1:nrow(est.states)))
  }
  
  if(is.sim){
    est.states$act.state.label <- het_cnt$state
    est.states$act.state <- NA
    est.states <- .genStateMap(est.states, est_cols = c('state', 'state.label'), 
                               act_cols = c('act.state', 'act.state.label'))
  }
  
  ret_obj <- est.states
  if(ret.raw){
    ret_obj <- list("df"=est.states, "model"=fit.mod)
    # metrics(ret_obj$df$state.label, ret_obj$df$act.state.label)$micro$F1
    # metrics(ret_obj$df$state.label, ret_obj$df$act.state.label)
    # summary(ret_obj$model)
    
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
  ret_obj <- list("map"=uniq_map, "model"=model)
  
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
    ret_obj[['act.map']] <- uniq_map_act
  }
  return(ret_obj)
}
