#' Multi-class performance metrics
#' @description Get multi-class performance metrics in the form of
#' Precision, Accuracy, and micro/macro F1 metrics. The base metrics
#' are calculated from the caret packages confusionMatrix() function
#' 
#' @param predictions vector of predictions
#' @param reference vector of actual classes
#'
#' @importFrom caret confusionMatrix
#' @importFrom assert_that assertthat
#' 
#' @return
#' 2 element list of a dataframe containing micro and macro F1 scores with 
#' precision and accuracy
#' @export
#'
#' @examples
#' metrics(predictions=sample(1:4, 1000, replace=T), reference=sample(1:4, 1000, replace=T))
metrics <- function(predictions, reference){
  # predictions <- hmm1$state
  # reference <- hmm1$act.state
  assert_that(length(predictions) == length(reference),
              msg="Predictions and reference are not the same vector length")
  
  # clean data
  levels <- unique(sort(c(predictions, reference)))
  predictions <- factor(predictions, levels=levels)
  reference   <- factor(reference, levels=levels)
  
  # CM
  cm <- confusionMatrix(predictions, reference)
  
  # microF1
  conf      <- as.data.frame(do.call(rbind, .confstats(cm)))
  pr        <- data.frame("Precision"=.pr(tp = sum(conf$tp), fp = sum(conf$fp)),
                          "Recall"=.re(tp = sum(conf$tp), fn = sum(conf$fn)))
  micro_f1  <- cbind(pr, "F1"=.f1(pr = pr$Precision, re = pr$Recall))
  
  # macroF1
  pr        <- cm$byClass[,c("Precision", "Recall")]
  f1        <- apply(pr, 1, function(i) .f1(pr = i['Precision'], re = i['Recall']))
  macro_f1  <- cbind(pr, "F1"=f1)
  
  return(list("cm"=cm$table, "micro"=micro_f1, "macro"=macro_f1))
}

.confstats <- function(cm){
  tbl <- cm$table
  
  conf <- lapply(colnames(tbl), function(ref_i){
    tp <- tbl[ref_i,ref_i]
    fp <- sum(tbl[ref_i, colnames(tbl) != ref_i])
    fn <- sum(tbl[colnames(tbl) != ref_i, ref_i])
    return(c("tp"=tp, "fp"=fp, "fn"=fn))
  })
  return(conf)
}

.pr <- function(tp, fp){
  tp / (tp + fp)
}

.re <- function(tp, fn){
  tp / (tp + fn)
}

.f1 <- function(pr, re){
  f1 <- 2 * ((pr * re) / (pr + re))
  if(is.nan(f1)) f1 <- 0
  f1
}
