#' Determine which associations are significant
#' 
#' Determines which associations are significant using emperical 
#' bayes and false discovery rate adjustment.
#' @param scores A symmetric p by p matrix of association scores.
#' @return A p by p matrix of likelihood.
#' @export
#' @title Determine which associations are significant
#' @description Determines which associations are significant using emperical 
#'   bayes and false discovery rate adjustment.
#' @param scores Either a vector or a symmetric matrix of association scores.
#' If a vector, `scores` should be the lower.tri entries of the the association
#' matrix.
#' @param gene_names (optional) Vector of gene names.
#' @param robust Should robust estimates of mean and variance be used?
#' @return A list containing the p by p matrix of likelihood ratios (fdr rates),
#' the estimate of mu for the null distribution, and the estimate of sigma
#' for the null distribution.
#' @export
#' @note For cPLS scores, first use `normalize_cpls_scores()` to normalize the
#' scores.
#' @note The robust estimators used are median and MAD.
fdr <- function(scores, gene_names = NULL, robust = TRUE) {
  
  # First obtain number of genes, `p`, and a vector of scores `val`.
  # Note: `scores` can be either a vector or matrix.
  if(length(dim(scores)) == 0) {
    # Vector of scores is provided.
    p <- sqrt(1 + 8*length(scores))/2 + 1/2 # Number of genes in network.
    if((p %% 1) != 0) stop("Length of scores does not match to lower.tri of matrix.")
    vals <- scores
  } else if(length(dim(scores)) == 2) {
    # Turn matrix of scores into vector.
    p <- nrow(scores) #Number of genes.
    vals <- scores[lower.tri(scores)] #Vector of scores between each gene.
    if(is.null(gene_names) & !is.null(colnames(scores))) {
      gene_names <- colnames(scores)
    }
  } else {
    stop("Scores should be a vector or 2 dimensional matrix.")
  }
  
  rm(scores)
  gc()
  
  # Estimate paramters of null distribution (assumed to be Gaussian).
  if(robust) {
    mu.f0 <- median(vals) # Use median since distribution of scores is skewed.
    sigma.f0 <- 1.4826 * median(abs(vals - median(vals))) # Use 1.4826 * MAD.
  } else {
    mu.f0 <- mean(vals) 
    sigma.f0 <- sd(vals) 
  }
  
  # Emperical Bayes FDR. Save likelihood ratios.
  likelihood <- matrix(0, nrow = p, ncol = p)
  likelihood[lower.tri(likelihood)] <- 
    dnorm(vals, mu.f0, sigma.f0) / approx(density(vals), xout = vals)$y
  likelihood <- likelihood + t(likelihood)
  diag(likelihood) <- Inf
  
  # Add gene names to likelihood matrix.
  colnames(likelihood) <- gene_names
  
  return(list(likelihood = likelihood, mu.f0 = mu.f0, sigma.f0 = sigma.f0))
}


#' Infer adjacency matrix from scores
#' 
#' Determines which associations are significant using emperical 
#' Bayes for false discovery rate adjustment.
#' @param scores A symmetric p by p matrix of association scores.
#' @param significance Must be specified if method = "fdr". 
#' Specifies the false discovery rate threshold to use.
#' @param method Character string specifing which method to use. Currently only
#' "fdr" is implemented.
#' @return The p by p adjacency matrix containing only 1's and 0's.
#' @export
get_adjacency_matrix <- function(scores, significance = 0.05, method = "fdr") {
  if(length(dim(scores)) != 2) {
    stop("scores should be a 2 dimensional matrix or data frame.")
  }
  if(is.null(method)) {
    stop("method must be specified.")
  }
  
  if(method == "fdr") {
    likelihood <- fdr(scores)$likelihood
    adj_matrix <- 1 * (likelihood < significance)
  }
  else {
    error("method must be one of: `fdr`.")
  }
  
  return(adj_matrix)
}