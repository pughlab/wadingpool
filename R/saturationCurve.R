#' Saturation curve estimations using asymptotic regression models
#' 
#' @description Calculates the nonlinear model fit using an 
#' Asymptotic Regression Model from and non-linear least squares
#' from the stats::nls and stats::SSasymp functions, where the
#' following parameters are being estimated:
#'  A = Response at 100% saturation (horizontal asymptote)
#'  R = Response when x = 0
#'  l = Natural logarithm of the rate constant
#'  
#'  Given an input of X values to predict on (pred), this function will 
#'  return the predicted Response based on the model. Or the 
#'  corresponding X for a given saturation (S), or a Saturation 
#'  associated with a given X.
#'
#' @param dat Two column data frame where dat[,1] is X values and dat[,2] is response
#' @param pred Series of X values to predict on (Default=NULL)
#' @param Xin X value to calculate the saturation point for (Default=NULL)
#' @param S Saturation point to calculate the X value for (Default=NULL)
#'
#' @import stats
#' @return
#' A list of the model data, coefficients, and returned values if pred, Xin, or
#' S are given
#' @export
saturationCurve <- function(dat, pred=NULL, S=NULL, Xin=NULL){
  # nlsdata <- data.frame("cov"=as.numeric(rownames(samp2)),
  #                       "SNPs"=samp2$n_mean)[-1,]
  # pred = seq(0,10,by=0.01)
  
  if(ncol(dat) != 2) stop("Input 'dat' must be a 2 column data.frame")
  
  colnames(dat) <- c("X", "Response")
  nlsfitSS <- nls(Response ~ SSasymp(X, Asym, R0, lrc),
                  data=dat)
  A <-summary(nlsfitSS)$coefficients[1]   # Response at 100% saturation (horizontal asymptote)
  R <-summary(nlsfitSS)$coefficients[2]   # Response when X = 0
  l <-summary(nlsfitSS)$coefficients[3]   # Natural logarithm of the rate constant
  
  ## Math Functions
  # Returns the X-value associated with a given saturation (S) value (S*A)
  getX <- function(A, R, l, S){ -1*( log(((S*A)-A)/(R-A)) / exp(l) ) }
  # Returns the saturation (S) associated with a given X value (Xin)
  getS <- function(A, R, l, Xin){ (exp((-1*Xin) * exp(l)) * (R-A) + A) / A }
  
  retobj <-list("model"=nlsfitSS,
                "Asymp_coef"=A,
                "R0_coef"=R,
                "lrc"=l)
  if(!is.null(pred)) retobj[['pred']] <- setNames(as.numeric(predict(nlsfitSS,list(X=pred))), pred)
  if(!is.null(S)) retobj[['Xfit']] <- setNames(getX(A, R, l, S), S)
  if(!is.null(S)) retobj[['Yfit']] <- setNames(as.numeric(predict(nlsfitSS,list(X=retobj[['Xfit']]))), S)
  if(!is.null(Xin)) retobj[['S']] <- setNames(getS(A, R, l, Xin), Xin)
  return(retobj)
}