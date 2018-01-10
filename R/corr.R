
#' Wrapper for correlation co-expression
#' 
#' Conducts co-expression analysis using correlation for association measure.
#' No preprocessing is done to x prior to computing the correlations.
#' @param x The n by p matrix of counts.
#' @param threshold Cutoff for significant associations. If NULL, all correlations
#' are returned. Otherwise, correlations of magnitude at or below this threshold are 
#' set to zero.
#' @param method The method used to compute correlations. Should be either "pearson"
#'  or "spearman". The default is Spearman, which provides 
#'  the more conservative estimation of associations.
#' @return A list containing the p by p matrix of correlations,
#' the threshold used to determine significant associations, and
#' the method used to compute the correlations.
#' @export
run_corr <- function(x, threshold = NULL, method = "spearman") {
  if(!(method %in% c("pearson", "spearman"))) {
    stop('method should be one of c("pearson", "spearman").')
  }
  
  if(is.integer(x[1, 1])) {
    x <- x + 0.0
  }
  scores <- cor(x, method = method)
  diag(scores) <- 0
  
  if(!is.null(threshold)) {
    scores[abs(scores) <= threshold] <- 0
  }

  return(list(scores = scores, threshold = threshold, method = method))
}
