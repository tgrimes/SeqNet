#' Wrapper for inverse covariance matrix
#' 
#' Conducts co-expression analysis using full partial correlations; these are
#' computed by inverting the sample correlation matrix of the observations.
#' @param x The n by p matrix of observations.
#' @param threshold Cutoff for significant associations. If NULL, all correlations
#' are returned. Otherwise, correlations of magnitude at or below this threshold are 
#' set to zero.
#' @param k An integer that ensures the matrix inverse is numerically stable. 
#' k = 1 is default; higher values will give less stable results.
#' @return A list containing the p by p matrix of partial correlations and
#' the threshold used to determine significant associations.
#' @export
run_icorr <- function(x, threshold = NULL, k = 1) {
  if(is.integer(x[1, 1])) {
    x <- x + 0.0
  }
  scores <- cor(x)
  scores[is.na(scores) | is.nan(scores)] <- 0
  
  # Add constant to diagonal to ensure PD, then solve.
  p <- nrow(scores)
  eigen_vals <- eigen(scores)$values
  scores <- scores + diag((max(eigen_vals) * 10^-k - min(eigen_vals)), p) * 
    (min(eigen_vals) < max(eigen_vals) * 10^-k) # Add diag if this is TRUE.
  
  scores <- corpcor::cor2pcor(scores) # Computes inverse and handles negative signs.
  
  if(!is.null(threshold)) {
    scores[abs(scores) <= threshold] <- 0
  }
  
  colnames(scores) <- colnames(x)
  
  return(list(scores = scores, 
              threshold = threshold))
}