#' Generate sample from overdispered Poisson
#' 
#' Generates count data based on the overdispered Poisson distribution.
#' @param n The number of samples to generate.
#' @param network The underlying network (from create_network())
#' @return A list containing the n by p matrix of samples, the underlying network,
#' and an n by p matrix containing mean expression values for each gene on each sample
#' after incorporating the underlying network.
#' @export
#' @examples
#' n <- 10
#' network <- create_network(p = 100, modules = list(1:100)) # Create small-world network.
#' mu <- get_reference_count_means()
#' # Generate RNA-seq data using the network and reference means:
#' x <- gen_gaussian(n, network, mu)$x
gen_gaussian <- function(n, 
                         network) {
  library(mvtnorm) #rmvnorm
  library(MCMCpack) #diwish() and riwish()
  
  if(n <= 0) stop("n must be positive.")
  p <- network$p

  graph <- get_adj_matrix_from_network(network)
  connected <- apply(graph, 2, function(x) sum(x != 0) > 0)
  

  x <- matrix(0, nrow = p, ncol = n) # Store samples in columns at first.
  A <- graph[connected, connected]
  n_A <- sum(lower.tri(A))
  A[lower.tri(A)] <- A[lower.tri(A)] * (-1)^rbinom(n_A, 1, 0.5) * runif(n_A, 0.5, 1)
  #A[lower.tri(A)] <- A[lower.tri(A)] * 0.75
  A[upper.tri(A)] <- 0
  A <- A + t(A)
  precision <- A + diag((1 + abs(min(eigen(A)$values))), nrow(A))
  eigen(precision)$values
  sigma <- solve(precision)
  x[connected, ] <- t(rmvnorm(n, sigma = sigma))
  x[!connected, ] <- rnorm(sum(!connected) * n)
  
  # x <- matrix(0, nrow = p, ncol = n) # Store samples in columns at first.
  # A <- graph[connected, connected]
  # n_A <- sum(lower.tri(A))
  # A[lower.tri(A)] <- A[lower.tri(A)] * (-1)^rbinom(n_A, 1, 0.5) * runif(n_A, 0.5, 1)
  # #A[lower.tri(A)] <- A[lower.tri(A)] * 0.75
  # A[upper.tri(A)] <- 0
  # A <- A + t(A)
  # precision <- A + diag((1 + abs(min(eigen(A)$values))), nrow(A))
  # eigen(precision)$values
  # sigma <- solve(precision)
  # 
  # for(i in 1:n) {
  #   x[connected, i] <- rmvnorm(1, sigma = sigma)
  #   x[!connected, i] <- rnorm(sum(!connected))
  # }
  
  x <- t(x) # Turn into n by p matrix.
  
  
  # Add column names to simulated dataset.
  if(is.null(network$node_names)) {
    colnames(x) <- paste(1:p)
  } else {
    colnames(x) <- network$node_names
  }
  
  if(any(is.na(x))) {
    warning(paste("NAs produced by `gen_gaussian`: replacing",
                  sum(is.na(x)), "of", length(x), "values with 0"))
    x[which(is.na(x))] <- 0
  }
  
  return(list(x = x, 
              network = network))
}

attr(gen_gaussian, "name") <- "gen_gaussian"


# # Check input for sigma.
# if(length(sigma) == 1) {
#   sigma <- rep(sigma, p)
# }
# if(length(sigma) != p) stop("sigma should be a vector of length 1 or p.")
# 
# # Check input for rho.
# if(length(rho) == 1) {
#   rho <- matrix(rho, p, p)
# }
# if(length(rho) != p^2) stop("rho should be a value or p by p matrix")
# diag(rho) <- 0

# @param sigma A nonegative value or vector of length p. Specifies the standard
# deviation in expression for each gene.
# @param rho A value or p by p matrix with entries between -1 and 1. Specifies
# the correlation between each gene. If two genes are not connected in the network,
# the correlation for that pair is ignored.
