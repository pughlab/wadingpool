#' Returns the tag from an execution
#' 
#' @description Looks for and parses the
#' 'Runtime: 25' line in shell executions
#'
#' @param results vector of printed statments
#' @param tag ID tag for the parsed output (e.g. Runtime)
#'
#' @return Data value for the parsed output
parseOutput <- function(results, tag){
  runtime <- results[grep(tag, results, ignore.case = TRUE)]
  runtime <- gsub("^.*: ", "", runtime)
  return(runtime)
}

#' rowxrow or colxcol execution
#'
#' @param mat Matrix 
#' @param fun Input function
#' @param margin Margin: col(2) or row(1)
matall <- function(mat, margin=2, fun){
  apply(mat, margin, function(i){
    apply(mat, margin, function(j){
      fun(i,j)
    })
  })
}
