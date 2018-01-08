
#' Wrapper for correlation co-expression
#' 
#' Conducts co-expression analysis using correlation for association measure.
#' @param x The n by p matrix of counts
#' @param threshold Cutoff for significant associations. If NULL, all correlations
#' are returned. Otherwise, correlations at or below this threshold are set to zero.
#' @return A list containing the p by p matrix of correlations,
#' and the threshold used to determine significant associations. 
#' @export
run_corr <- function(x, threshold = 0.9) {
  if(is.integer(x[1, 1])) {
    x <- x + 0.0
  }
  scores <- cor(x)
  diag(scores) <- 0
  
  if(!is.null(threshold)) {
    scores[scores <= threshold] <- 0
  }

  return(list(scores = scores, threshold = threshold))
}
