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
  switch.val <- 100 * round(switch.prob, 2)
  rand.val <- sample(1:100, N, replace = T)
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
#' @param x Numeric vector
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
sanitizeData <- function(x){
  x <- log(x / median(x))
  x <- Winsorize(x, probs=c(0.01, 0.99))
  return(x)
}