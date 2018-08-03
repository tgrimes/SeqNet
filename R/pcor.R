
#' Wrapper for partial correlations from corpcor
#' 
#' Conducts co-expression analysis using partial correlations for association 
#' measure.
#' @param x The n by p matrix of counts.
#' @param threshold Cutoff for significant associations. If NULL, all correlations
#' are returned. Otherwise, correlations of magnitude at or below this threshold are 
#' set to zero.
#' @return A list containing the p by p matrix of partial correlations,
#' the threshold used to determine significant associations, and
#' the method used to compute the correlations.
#' @export
run_pcor <- function(x, threshold = NULL, verbose = FALSE) {
  if(is.integer(x[1, 1])) {
    x <- x + 0.0
  }
  
  scores <- pcor.shrink(x, verbose = verbose)
  
  if(!is.null(threshold)) {
    scores[abs(scores) <= threshold] <- 0
  }
  
  colnames(scores) <- colnames(x)
  
  return(list(scores = scores, threshold = threshold))
}
