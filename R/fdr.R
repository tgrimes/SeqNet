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
#' @param threshold The threshold for false discovery rate, used to determine 
#' significance; scores with ikelihood ratios above this threshold are
#' set to zero. If NULL, all scores are preserved.
#' @param transformation (option) Function applied to scores before
#' estimating the null distribution. If provided, the function should take 
#' a vector as input and return a vector for output.
#' @param robust Should robust estimates of mean and variance be used?
#' @param gene_names (optional) Vector of gene names.
#' @param include_likelihood Should the matrix of likelihood ratios be provided
#' in the output?
#' @return A list containing the scores with non-significant values set to zero,
#' the estimate of mu for the null distribution, the estimate of sigma
#' for the null distribution, the empirical density estimate, and 
#' (optional) likelihood ratios (fdr rates) for each score.
#' @export
#' @note For cPLS scores, it is recommended to use `normalize_cpls_scores()`
#' as the transformation.
#' @note The robust estimators used are median and MAD.
fdr <- function(scores, 
                threshold = 0.05,
                transformation = NULL, 
                robust = TRUE,
                include_likelihood = FALSE) {
  
  # Determine whether scores is a vector or matrix.
  if(is.matrix(scores)) {
    if(length(dim(scores)) != 2) stop("scores is not a 2 dimensional matrix.")
    if(!all(scores == t(scores))) stop("scores is not symmetric.")
    score_names <- colnames(scores)     # Keep column names for scores.
    scores <- scores[lower.tri(scores)] # Turn scores into vector.
    matrix_provided <- TRUE
  } else if(is.vector(scores)) {
    matrix_provided <- FALSE
  } else {
    stop("scores is neither a vector nor symmetric matrix.")
  }
  
  # Apply transformation, if one is provided.
  if(!is.null(transformation)) {
    scores <- transformation(scores)
    if(!is.vector(scores)) stop("transformation did not return a vector as output.")
  }
  
  # Estimate paramters of null distribution (assumed to be Gaussian).
  if(robust) {
    mu.f0 <- median(scores) # Use median to account for any skewness.
    sigma.f0 <- 1.4826 * median(abs(scores - median(scores))) # Use 1.4826 * MAD.
  } else {
    mu.f0 <- mean(scores) 
    sigma.f0 <- sd(scores) 
  }
  
  # Emperical Bayes FDR. Save likelihood ratios.
  f <- density(scores)
  likelihood <- dnorm(scores, mu.f0, sigma.f0) / approx(f, xout = scores)$y
  
  # Apply thresholding if provided.
  if(!is.null(threshold)) {
    scores[likelihood > threshold] <- 0
  }
  
  if(matrix_provided) {
    # Put scores back into matrix form.
    scores <- vector_to_matrix(scores, diag_val = 0)
    colnames(scores) <- score_names
    
    # If likelihoods are requested, put them into matrix form.
    if(include_likelihood) {
      likelihood <- vector_to_matrix(likelihood, diag_val = Inf)
      colnames(likelihood) <- score_names
    }
  } else {
    if(include_likelihood) {
      names(likelihood) <- names(scores) # Pair score names with likelihoods.
    }
  }
  
  fdr_return_list <- list(scores = scores, 
                          mu.f0 = mu.f0, 
                          sigma.f0 = sigma.f0,
                          f = f)
  if(include_likelihood) {
    fdr_return_list <- c(fdr_return_list, 
                         list(likelihood = likelihood))
  }
  
  return(fdr_return_list)
}


test_fdr <- function() {
  # Note: there is additional code throughout for testing the performance of
  # the fdr procedure. This can be expanded on in the future.
  
  n_null <- 1000
  n_alt <- 100
  scores <- rnorm(n_null, 0, 1)
  scores <- c(scores, rnorm(n_alt / 2, 4, 1), rnorm(n_alt / 2, -4, 1))
  index_null <- 1:n_null
  index_alt <- (n_null + 1):(n_null + n_alt)
  
  # Test with vector input.
  s <- fdr(scores, robust = FALSE) # Check if robust = FALSE throws error.
  s <- fdr(scores)
  cat("Vector input returns vector output:", is.vector(s$scores), "\n")
  mean(s$scores[index_null] == 0) # specificity
  mean(s$scores[index_alt] != 0) # sensitivity
  1 - sum(which(s$scores != 0) %in% index_alt) / sum(s$scores != 0) # fdr
  
  # Test with matrix input.
  scores <- matrix(rnorm(100^2), 100, 100)
  scores <- scores + t(scores)
  s <- fdr(scores, robust = FALSE)
  s <- fdr(scores)
  cat("Matrix input returns matrix output:", is.matrix(s$scores), "\n")
  sum(s$scores != 0) # Number of significant scores.

  # Test using transformation.
  n_null <- 1000
  n_alt <- 100
  scores <- rlnorm(n_null, 0, 1)
  scores <- c(scores, rlnorm(n_alt, 4, 1))
  index_null <- 1:n_null
  index_alt <- (n_null + 1):(n_null + n_alt)
  s <- fdr(scores, transformation = log)
  mean(s$scores[index_null] == 0) # specificity
  mean(s$scores[index_alt] != 0) # sensitivity
  1 - sum(which(s$scores != 0) %in% index_alt) / sum(s$scores != 0) # fdr
  
  s <- fdr(scores, include_likelihood = TRUE)
  cat("Likelihood is returned when requested:", !is.null(s$likelihood), "\n")
  
  s <- fdr(scores, threshold = NULL)
  cat("NULL threshold preserves all scores:", !any(s$scores == 0), "\n")
}

#' Turn a vector into a matrix
#' 
#' The input vector is assumed to be the lower.tri of a symmetric matrix. This
#' function recreates the original matrix from the vector.
#' @param x the vector containing lower.tri entries of a symmetric matrix.
#' @param diag_val value to put along the diagonal.
#' @return A symmetric matrix.
#' @export
vector_to_matrix <- function(x, diag_val = 0) {
  if(!is.vector(x)) stop("x should be a vector.")
  
  p <- sqrt(1 + 8 * length(x)) / 2 + 0.5 # Dimension of original matrix.
  if((p %% 1) != 0) stop("x is not the lower.tri entries of a square matrix.")
  
  y <- matrix(0, p, p) # Setup the output matrix
  y[lower.tri(y)] <- x # Add x to the lower.tri
  y <- y + t(y)        # Add x to the upper.tri (maintaining symmetry).
  
  diag(y) <- diag_val
  
  return(y)
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