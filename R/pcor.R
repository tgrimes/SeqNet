#' Wrapper for partial correlations from corpcor
#' 
#' Conducts co-expression analysis using full partial correlations; these are
#' computed using the shrinkage approach for covariance estimation from the 
#' corpcor package.
#' @param x The n by p matrix of counts.
#' @param threshold Cutoff for significant associations. If NULL, all correlations
#' are returned. Otherwise, correlations of magnitude at or below this threshold are 
#' set to zero.
#' @param verbose If true, progress messages will be printed from pcor.shrink().
#' @return A list containing the p by p matrix of partial correlations and
#' the threshold used to determine significant associations.
#' @export
run_pcor <- function(x, threshold = NULL, verbose = FALSE) {
  if(is.integer(x[1, 1])) {
    x <- x + 0.0
  }
  
  scores <- corpcor::pcor.shrink(x, verbose = verbose)
  
  if(!is.null(threshold)) {
    scores[abs(scores) <= threshold] <- 0
  }
  
  colnames(scores) <- colnames(x)
  
  return(list(scores = scores, 
              threshold = threshold))
}


#' Test if a shrinkage object is symmetric.
#' 
#' @param x An object of class 'shrinkage'.
#' @param ... Further argumnts passed to [isSymmetric()].
#' @return Logical indicated if x is symmetric or not.
#' @export
isSymmetric.shrinkage <- function(x, ...) {
  class(x) <- "matrix"
  return(isSymmetric(x, ...))
}

