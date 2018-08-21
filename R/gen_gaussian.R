#' Generate sample from a Gaussian graphical model.
#' 
#' Generates count data based on the multivariate normal distribution.
#' @param n The number of samples to generate.
#' @param network The underlying network (from [create_network()]). 
#' @param sigma The covariance matrix to use. If network$sigma or network$precision 
#' is defined, they will be used instead of sigma. Otherwise, if sigma NULL, 
#' a matrix is generated for sigma.
#' @param network_list A list of networks (from [create_network()]).
#' @return A list containing the n by p matrix of samples, the underlying network,
#' and the p by p covariance matrix.
#' @export
#' @examples
#' n <- 10
#' # Create a simple small-world network.
#' network <- create_network(p = 100, modules = list(1:100)) 
#' # Generate "normalized" expression data from the network:
#' x <- gen_gaussian(n, network)$x
gen_gaussian <- function(n, network = NULL, sigma = NULL, network_list = NULL) {
  if(!is.null(network_list)) 
    return(gen_gaussian_list(n = n, network_list = network_list))
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
      sigma <- get_sigma_from_precision_list(list(precision))[[1]]
    } else {
      precision <- random_precision_from_network_list(list(network))[[1]]
      sigma <- get_sigma_from_precision_list(list(precision))[[1]]
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


#' Generate samples from a Gaussian graphical model.
#' 
#' Generates count data based on the multivariate normal distribution. Input
#' consists of a list of networks, for which corresponding samples are generated.
#' @param n The number of samples to generate.
#' @param network_list A list of networks (from [create_network()]).
#' @return A list containing output from [gen_gaussian()] for each network.
#' @export
gen_gaussian_list <- function(n, network_list) {
  precision_list <- random_precision_from_network_list(network_list)
  sigma_list <- get_sigma_from_precision_list(precision_list)
  gen_list <- lapply(sigma_list, 
                     function(sigma) gen_gaussian(n = n, sigma = sigma))
  for(i in 1:length(gen_list)) {
    gen_list[[i]]$network <- network_list[[i]]
  }
  return(gen_list)
}



#' Generate partial correlations for a list of networks.
#' 
#' Matrices containing partial correlations are generated based on the structure 
#' of each network. Edges that are common across networks are given the same 
#' partial correlation value.
#' @param network_list A list of network objects (from [create_network()]).
#' @param k An integer that ensures the matrix inverse is numerically stable. 
#' k = 1 is default; higher values will give less stable results.
#' @param limits A vector of length 2 containing the lower and upper limits
#' for partial correlations (generated from a uniform distribution).
#' @return A list of matrices containing partial correlations.
#' @export
random_precision_from_network_list <- function(network_list, k = 1, 
                                               limits = c(0.5, 1)) {
  if(length(limits) != 2) 
    stop("limits must be a vector of length 2.")
  if(limits[1] > 1 | limits[1] < 0)
    stop("limits[1] is not between 0 and 1.")
  if(limits[2] > 1 | limits[2] < 0)
    stop("limits[2] is not between 0 and 1.")
  if(limits[1] >= limits[2])
    stop("limits[1] must be less than limits[2]")
  
  nodes <- unique(unlist(lapply(network_list, function(nw) nw$node_names)))
  p <- length(nodes)
  
  # Generate partial correlation values for every possible connection.
  precision_all <- matrix(0, p, p)
  m <- p * (p - 1) / 2
  precision_all[lower.tri(precision_all)] <- 
    (-1)^rbinom(m, 1, 0.5) * runif(m, limits[1], limits[2])
  precision_all <- precision_all + t(precision_all)
  
  # Obtain an adjacency matrix representation of each network.
  adj_matrix_list <- lapply(network_list, get_adj_matrix_from_network)
  
  # Fill in partial correlations for each network.
  precision_list <- lapply(adj_matrix_list, function(nw) {
    index <- which(nodes %in% colnames(nw))
    
    # Create p by p precision matrix with partial correlations.
    precision <- nw * precision_all[index, index]
    diag(precision) <- 1
    return(precision)
  })
  
  # Ensure each precision matrix is invertible by adjusting the diagonal.
  eigen_val_list <- lapply(precision_list, function(P) eigen(P)$values)
  adjustment_list <- lapply(eigen_val_list, function(E) {
    (max(E) * 10^-k - min(E)) * (min(E) < max(E) * 10^-k)
  })
  # Apply the same adjustment to each precision matrix.
  # This ensures the partial correlations remain constant for each connection.
  adjustment <- diag(max(unlist(adjustment_list)), p)
  precision_list <- lapply(precision_list, function(P) P + adjustment)
  
  # Standardize with 1s along diagonal so that partial correlation = -Precision.
  precision_list <- lapply(precision_list, cov2cor)
  
  return(precision_list)
}


#' Compute covariance matrices from a list of precision matrices
#' 
#' @param precision_list A list of precision matricies.
#' @return A list of covariance matrices. 
#' @export
get_sigma_from_precision_list <- function(precision_list) {
  # If a network is given, check if "precision" matrix is provided.
  if(class(precision_list[[1]]) == "network") {
    if(!all(sapply(precision_list, function(nw) "precision" %in% names(nw))))
      stop("List of network objects provided but not all contain precision matrix.")
    precision_list <- lapply(precision_list, function(nw) nw$precision)
  }
  
  if(!all(sapply(precision_list, is.matrix))) 
    stop("At least one object in precision list is not a matrix.")
  if(!all(sapply(precision_list, function(P) isSymmetric(unname(P))))) 
    stop("At least one matrix in precision list is not symmetric.")
  
  # Invert adjusted precision matrix to obtain covariance matrix.
  sigma_list <- lapply(precision_list, function(P) pcor2cor(P))
  
  return(sigma_list)
}
