
#' Wrapper for inverse correlation matrix
#' 
#' Conducts co-expression analysis using partial correlations for association 
#' measure, computed by inverting the correlation matrix.
#' @param x The n by p matrix of counts.
#' @param threshold Cutoff for significant associations. If NULL, all correlations
#' are returned. Otherwise, correlations of magnitude at or below this threshold are 
#' set to zero.
#' @param method The method used to compute correlations. Should be either "pearson"
#'  or "spearman". The default is Spearman, which provides 
#'  the more conservative estimation of associations.
#' @return A list containing the p by p matrix of partial correlations,
#' the threshold used to determine significant associations, and
#' the method used to compute the correlations.
#' @export
run_icorr <- function(x, threshold = NULL, method = "pearson") {
  if(!(method %in% c("pearson", "spearman"))) {
    stop('method should be one of c("pearson", "spearman").')
  }
  
  if(is.integer(x[1, 1])) {
    x <- x + 0.0
  }
  scores <- cor(x, method = method)
  scores[is.na(scores) | is.nan(scores)] <- 0
  scores <- scores + diag((1 + abs(min(eigen(scores)$values))), nrow(scores))
  scores <- solve(scores)
  
  if(!is.null(threshold)) {
    scores[abs(scores) <= threshold] <- 0
  }
  
  colnames(scores) <- colnames(x)
  
  return(list(scores = scores, threshold = threshold, method = method))
}
