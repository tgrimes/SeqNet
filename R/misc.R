#' Check if a matrix is positive definite
#' 
#' @param x A matrix to check.
#' @return Returns TRUE if the matrix is positive definite and FALSE otherwise.
#' @export
is_PD <- function(x) {
  # Note: for matricies with p > 500, computing the smallest eigenvalue is 
  # faster than Choleski. See RSpectra::eigs_sym(adj, k = 1, which = "SA")
  
  # Check diagonal of Choleski factorization. If the matrix is not PD, the
  # test will fail (i.e. be FALSE).
  fails_test <- tryCatch(any(diag(chol(x)) < 0), 
                         error = function(e) TRUE)
  return(!fails_test)
}