#' Set low ranking scores to zero.
#' 
#' The association scores are ranked (in magnitude) from strongest to weakest, 
#' and the bottom ranking scores are set to zero. This is useful when plotting,
#' if it is desired to show only the edges with the strongest associations.
#' @param scores The matrix of association scores.
#' @param rank Scores below this rank are set to zero.
#' @param verbose If true, the number of associations that were set to zero 
#' are printed.
#' @return The matrix of association scores with low ranking scores set to zero. 
#' Note, in the case of ties (many edges having the same association), the total
#' number of edges may exceed rank.
#' @export
set_low_rank_scores_to_zero <- function(scores, rank = 100, verbose = FALSE) {
  
  # Find descending order (in absolute value).
  index_top <- order(abs(scores), decreasing = TRUE)
  
  # Find the value at the cutoff index. 
  # Note 1: The value is used, rather than the index itself, so that ties can
  #         be handled appropriately.
  # Note 2: Since scores is symmetric, each value is duplicated, so twice of 
  #         rank is used.
  threshold <- abs(scores[index_top[2 * rank]])
  
  # Determine which scores are below the threshold.
  index_low_scores <- which(abs(scores) < threshold)
  
  if(verbose) {
    cat("Setting scores below rank", rank, "to zero.\n")
    num_ties <- sum(scores == threshold) / 2
    if(num_ties != 1) {
      cat("\tThere are", num_ties, "genes tied for rank", rank, ".\n")
    }
    cat("\t", sum(scores[index_low_scores] != 0) / 2, "of", 
        length(index_low_scores) / 2, "scores were non-zero.\n")
  }
  
  scores[index_low_scores] <- 0 # Set scores below top rank to 0. 
  return(scores)
}