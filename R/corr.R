
#' Wrapper for correlation co-expression
#' 
#' Conducts co-expression analysis using correlation for association measure.
#' @param x The n by p matrix of counts
#' @param threshold Cutoff for significant associations
#' @return A list containing `scores`, a p by p matrix of correlations, and 
#' `adj_matrix`, a p by p adjacency matrix.
#' @export
run_corr <- function(x, threshold = 0.9) {
  if(is.integer(x[1, 1])) {
    x <- x + 0.0
  }
  correlations <- cor(x)
  diag(correlations) <- 0
  
  x[which(is.na(x))] <- 0
  
  adj_matrix <- 1 * (correlations > threshold)
  
  return(list(scores = correlations, adj_matrix = adj_matrix))
}
