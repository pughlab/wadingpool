#' simulate
#' @description Simulates HMM data using a given family of distributions, 
#' an input number of states, and a transition probability.
#' 
#' @param N Integer value for number of elements to randomly generate
#' @param cnts Named integer vector for mean/lambda value for that state
#' @param switch.prob Numeric value for transition probability (Default=0.1)
#' @param family Family of distributions (Default=poission)
#'
#' @return
#' A Dataframe consisting of states and the simulated observations from each state
#' @export
#'
#' @examples
#' set.seed(1234)
#' simulate(100, cnts=setNames(c(10, 20), c('loh', 'het')), 
#'          switch.prob=0.05, family='poisson')
simulate <- function(N, cnts=c('loh'=10, 'het'=15), switch.prob = 0.10, family='poisson'){
  switch.val <- 10000 * round(switch.prob, 5)
  rand.val <- sample(1:10000, N, replace = T)
  sim_cnts <- switch(family,
                     poisson=lapply(cnts, function(cnt) rpois(N,cnt)),
                     gaussian=lapply(cnts, function(cnt) rnorm(n=N, mean=cnt, sd=sqrt(cnt))))
  stopifnot(!is.null(names(cnts)))
  all_states <- names(cnts)
  
  # states 
  draws <- data.frame(state = rep(NA, N), 
                      obs = rep(NA, N), 
                      dice = rand.val)
  
  # Set initial state
  init_state <- sample(all_states,1)
  draws$state[1] <- init_state
  draws$obs <- sim_cnts[[init_state]][1]
  
  for(k in 2:N){
    state_k <- draws$state[k-1]
    if(draws$dice[k-1] <= switch.val){
      # switch to a random state
      state_k <- sample(all_states[grep(state_k, all_states, invert = T)], 1)   # random state
    }
    draws$state[k] <- state_k
    draws$obs[k] <- sim_cnts[[state_k]][k]
  }
  # return
  return(cbind(roll = 1:N, draws))
}

#' Sanitize the data
#' @description If data is distributed according to a poisson distribution, median normalize
#' the data and log transform to make it more gaussian like
#'
#' @param x Numeric vector
#' @param winsorize_p Numeric value for probabilty to winsorize at (Default = 0.01)
#' 
#' @importFrom DescTools Winsorize
#' @importFrom stats median
#'
#' @return
#' Numeric vector of median normalized x in log space, followed by a 0.01,0.099 winsorization 
#' @export
#'
#' @examples
#' sanitizeData(rpois(100, 10))
sanitizeData <- function(x, winsorize_p=0.01){
  x <- log(x / median(x))
  
  ## Find minimum quantile that does not return an infinite value
  lowq <- 0; highq <- 1
  while(is.infinite(quantile(x, lowq))){
    lowq <- lowq + 0.01
  }
  while(is.infinite(quantile(x, highq))){
    highq <- highq - 0.01
  }
  q <- max(c(winsorize_p, lowq, 1-highq))
  if(q > winsorize_p){
    warning(paste0("Winsorize value of ", winsorize_p, " returns an infinite value, defaulting \n",
                 "to ", q, " which will return the next best winsorization probability"))
  }
  
  x <- Winsorize(x, probs=c(q, 1-q))
  return(x)
}


#' Parse CN seg files
#' @description Parses a seg file with given chromosome-separated
#' intervals and a discrete CN state. It will identify the number of
#' states occupying `state_frac` of the genome, where each interval
#' is of at least `min_frac` in size.
#' 
#' @param segf path to CN seg file
#' @param segcol character: column name for discrete CN state
#' @param start character: column name for interval start
#' @param end character: column name for interval end
#' @param state_frac numeric: minimal fraction of the genome occupied by CN states (default=0.90)
#' @param min_frac numeric: minimum segment size in fractions (default=0.05)
#' @param get_tt boolean: return a transition state probability matrix (in.dev)
#' @param get_seg boolean: return the sef file parsed
#'
#' @importFrom assertthat assert_that
#' @return
#' A one or two element list:
#'   'state': integer number of states
#'   'tt': transition matrix
#' @export
parseSeg <- function(segf, segcol='event', start='start', end='end', 
                     state_frac=0.9, min_frac=0.05, get_tt=FALSE, get_seg=FALSE){
  # segf <- 'net-001a.seg'
  assert_that(file.exists(segf), msg="Seg file does not exist")
  seg <- read.table(segf, header = TRUE, sep="\t", stringsAsFactors = FALSE, check.names = FALSE)
  # .chk_seg(seg, st_end=c(start, end), event=segcol)    # to be implemented
  
  # Get genomic size of all events
  seg$width   <- seg[,end] - seg[,start]
  segl_width  <- split(seg$width, seg[,segcol])
  seg_frac    <- sapply(segl_width, sum)/sum(seg$width)
  
  # Identify number of states that occupy f=`state_frac` of the genome
  small_segs  <- which(seg_frac < min_frac)
  max_frac    <- sum(seg_frac[-small_segs])  # max for all states > `min_frac`
  cumseg      <- cumsum(sort(seg_frac[-small_segs], decreasing = T))      # cumu
  state_cut   <- which(diff(cumseg > state_frac) == 1) + 1
  
  # Estimate transition prob
  ## tileGenome() to split the genome in the same way as hetCnts
  ## findOverlap() to overlap seg and tiles
  ## Generate expected transition matrix 
  ### tt <- table(x[-length(x)], x[-1])
  ### tt <- tt / rowSums(tt)
  
  ret_obj <- list("states"=as.integer(state_cut))
  if(get_tt){
    stop(paste0("Obtaining the transition matrix has not been implemented yet"))
    # ret_obj[['tt']] <- tt
  } else if(get_seg){
    ret_obj[['seg']] <- seg
  }
  return(ret_obj)
}

#' GRanges to Bed
#' @description simple function to convert a GRanges object to a data
#' frame in .bed file format
#' @param gr GRanges object
#' @importFrom GenomeInfoDb seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges strand
#' 
#' @return
#' Dataframe
#' @export
toBed <- function(gr){
  df <- data.frame(chrom=seqnames(gr),
                   chromStart=as.integer(start(gr)),
                   chromEnd=as.integer(end(gr)),
                   name=paste0("bin", c(1:length(gr))),
                   score=c(rep(".", length(gr))),
                   strand=strand(gr))
  return(df)
}


#' Zygosity combiner helper
#' @description Combines the zygosity data with the zygosity
#' position
#' @param zyg data.frame object returned from fitHMM()$df
#' @param zygpos GRanges object that corresponds with zyg
#' @param ret return 'reduced' or a 'raw' state (bins or segments)
#' @importFrom S4Vectors mcols<-
#' @importFrom S4Vectors mcols
#' @importFrom IRanges reduce
#' @return 
#' A GRanges object that maps zygpos to zyg
#' @export
combZygPos <- function(zyg, zygpos, ret='reduced'){
  # zyg <- hmm1$df
  # zygpos <- tiles_ov
  
  # Uncollapsed data
  seg <- zygpos[zyg$bin]
  mcols(seg) <- zyg
  print(head(seg))
  
  # Collapsed seg states
  if(ret=='reduced'){
    seg <- unlist(IRanges::reduce(split(seg, seg$state.label)))
    seg <- sort(seg)
    seg$state <- names(seg)
  }
  return(seg)
}