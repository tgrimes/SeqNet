#' Generate sample from a Gaussian graphical model.
#' 
#' Generates count data based on the multivariate normal distribution.
#' @param n The number of samples to generate.
#' @param network The underlying network (from [create_network()]). 
#' @param sigma The covariance matrix to use. If network$sigma or network$precision 
#' is defined, they will be used instead of sigma. Otherwise, if sigma NULL, 
#' a matrix is generated for sigma.
#' @return A list containing the n by p matrix of samples, the underlying network,
#' and the p by p covariance matrix.
#' @export
#' @examples
#' n <- 10
#' # Create a simple small-world network.
#' network <- create_network(p = 100, modules = list(1:100)) 
#' # Generate "normalized" expression data from the network:
#' x <- gen_gaussian(n, network)$x
gen_gaussian <- function(n, network = NULL, sigma = NULL) {
  if(is.null(network) & is.null(sigma))
    stop("Either network or sigma must be provided.")
  if(n <= 0) 
    stop("n must be positive.")
  
  # If network is provided, create covariance matrix sigma.
  if(!is.null(network)) {
    if(("sigma" %in% names(network)) & !is.null(network$sigma)) {
      sigma <- network$sigma
    } else if(("precision" %in% names(network)) & !is.null(network$precision)) {
      precision <- network$precision
      sigma <- get_sigma_from_precision(precision)
    } else {
      precision <- random_precision_from_network(network)
      sigma <- get_sigma_from_precision(precision)
    }
    gene_names <- network$node_names
  } else {
    if(ncol(sigma) != nrow(sigma))
      stop("sigma is not a square matrix.")
    if(any(diag(sigma) < 0))
      stop("sigma contains negative values along its diagonal.")
    gene_names <- colnames(sigma)
  }
  
  # Generate samples.
  x <- mvtnorm::rmvnorm(n, sigma = sigma)
  
  # Add column names to simulated dataset.
  if(is.null(gene_names)) {
    colnames(x) <- paste(1:ncol(x))
  } else {
    colnames(x) <- gene_names
  }
  
  if(any(is.na(x))) {
    warning(paste("NAs produced by `gen_gaussian`: replacing",
                  sum(is.na(x)), "of", length(x), "values with 0"))
    x[which(is.na(x))] <- 0
  }
  
  return(list(x = x, 
              sigma = sigma,
              network = network))
}

attr(gen_gaussian, "name") <- "gen_gaussian"


#' Generate a precision matrix from a network
#' 
#' A random precision matrix (i.e. an inverse covariance matrix) is generated
#' based on the structure from a given network.
#' @param network The underlying network (from [create_network()]). 
#' @return A p by p precision matrix.
#' @export
random_precision_from_network <- function(network) {
  # Obtain an adjacency matrix representation of the network.
  graph <- get_adj_matrix_from_network(network)
  p <- ncol(graph)
  
  # Find nodes with edges in the graph.
  connected <- apply(graph, 2, function(x) sum(x != 0) > 0)
  n_connected <- sum(connected)
  
  # Each edge corresponds to a nonzero partial correlation.
  # Generate random value for each partial correlation.
  A <- graph[connected, connected] # Subset graph on nodes with edges.
  n_A <- n_connected * (n_connected - 1) / 2 # Number of possible edges in the graph.
  A[lower.tri(A)] <- A[lower.tri(A)] * (-1)^rbinom(n_A, 1, 0.5) * runif(n_A, 0.75, 1)
  A[upper.tri(A)] <- 0
  A <- A + t(A)
  
  # Create p by p precision matrix with partial correlations.
  precision <- matrix(0, p, p)
  precision[connected, connected] <- A
  diag(precision) <- 1
  
  return(precision)
}

#' Compute the covariance matrix from a precision matrix
#' 
#' If the precision matrix is singular, it is made invertible by adding a 
#' diagonal matrix.
#' @param precision The precision matrix.
#' @param k An integer that ensures the matrix inverse is numerically stable. 
#' k = 1 is default; higher values will give less stable results.
#' @return A p by p covariance matrix.
#' @export
get_sigma_from_precision <- function(precision, k = 1) {
  # If a network is given, check if "precision" matrix is provided.
  if(is.list(precision)) {
    if(!("precision" %in% names(precision)))
      stop("precision matrix not found.")
    precision <- precision$precision
  }
  
  if(!is.matrix(precision)) stop("precision is not a matrix.")
  if(!isSymmetric(precision)) stop("precision matrix is not symmetric.")
  
  # Add values to diagonal to ensure positive definiteness and 
  # numerical stablility of taking the inverse.
  p <- nrow(precision)
  eigen_val <- eigen(precision)$values
  precision <- precision + diag((max(eigen_val) * 10^-k - min(eigen_val)), p) * 
    (min(eigen_val) < max(eigen_val) * 10^-k)
  
  # Invert precision matrix to obtain covariance matrix.
  sigma <- solve(precision)
  
  return(sigma)
}
