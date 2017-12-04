#' Determine which associations are significant
#' 
#' Determines which associations are significant using emperical 
#' bayes and false discovery rate adjustment.
#' @param scores A symmetric p by p matrix of association scores.
#' @return A p by p matrix of likelihood.
#' @export
fdr <- function(scores) {
  if(length(dim(scores)) != 2) {
    stop("scores should be a 2 dimensional matrix or data frame.")
  }
  
  p <- nrow(scores) #Number of genes.
  vals <- scores[lower.tri(scores)] #Vector of scores between each gene.
  
  mu.f0 <- mean(vals)
  sigma.f0 <- sd(vals)
  
  #Emperical Bayes FDR. Save likelihood of each observation.
  likelihood <- matrix(0, nrow = p, ncol = p)
  likelihood[lower.tri(likelihood)] <- 
    dnorm(vals, mu.f0, sigma.f0) / approx(density(vals), xout = vals)$y
  likelihood <- likelihood + t(likelihood)
  diag(likelihood) <- Inf
  
  return(list(likelihood = likelihood))
}

#' Infer adjacency matrix from scores
#' 
#' Determines which associations are significant using emperical 
#' Bayes for false discovery rate adjustment.
#' @param scores A symmetric p by p matrix of association scores.
#' @param significance False discovery rate to use.
#' @return A list containing the p by p adjacency matrix and all significant 
#'   scores.
#' @export
get_adjacency_matrix <- function(scores, significance = 0.1, method = 1) {
  if(length(dim(scores)) != 2) {
    stop("scores should be a 2 dimensional matrix or data frame.")
  }
  
  p <- nrow(scores) #Number of genes.
  vals <- scores[lower.tri(scores)] #Vector of scores between each gene.
  
  mu.f0 <- mean(vals)
  sigma.f0 <- sd(vals)
  
  ordered_vals <- order(vals)
  n_vals <- length(vals)
  
  #Emperical Bayes FDR.
  numerator <- dnorm(vals, mu.f0, sigma.f0)
  denominator <- approx(density(vals), xout = vals)$y
  plot(vals[ordered_vals], numerator[ordered_vals], col = "blue", type = "l",
       main = "Empirical Bayes fdr")
  lines(vals[ordered_vals], denominator[ordered_vals], col = "orange")
  lines(vals[ordered_vals], numerator[ordered_vals] / significance, col = "red")
  likelihood <- numerator / denominator
  
  if(method == 1) significant_scores <- likelihood < significance
  
  #Find significant edges.
  adj_matrix <- matrix(0, nrow = p, ncol = p)
  adj_matrix[lower.tri(adj_matrix)] <- 1 * significant_scores
  adj_matrix <- adj_matrix + t(adj_matrix)
  
  return(adj_matrix)
}