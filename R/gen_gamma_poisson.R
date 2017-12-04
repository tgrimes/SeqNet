# Gamma-Poisson model

# Requires: 
#   SimSeq - for default reference dataset
# library(SimSeq)

#' Obtain empirical distribution of mean RNA-seq counts
#'
#' The average expression level of each gene in a reference dataset are 
#' used as an empirical distribution. If no reference dataset is provided,
#' the kidney data from the SimSeq package are used.
#' @param reference_dataset A reference dataset of RNA-seq expression counts. 
#' Can be a matrix or data.frame object.
#' @param min_mean_count Average counts below this value will be excluded.
#' @param max_mean_count Average counts above this value will be excluded.
#' @return A vector containing the average count for each gene.
#' @export
#' @examples get_reference_means()
get_reference_count_means <- function(reference_dataset = NULL,
                                      min_mean_count = 30,
                                      max_mean_count = NULL) {
  if(is.null(reference_dataset)) {
    data(kidney)
    reference_data <- kidney$counts
  }
  avg_counts <- apply(reference_data, 1, mean)
  
  if(!is.null(min_mean_count)) {
    avg_counts <- avg_counts[avg_counts >= min_mean_count]
  }
  
  if(!is.null(max_mean_count)) {
    avg_counts <- avg_counts[avg_counts <= max_mean_count]
  }
  
  return(avg_counts)
}


#' Generate sample from overdispered Poisson
#' 
#' Generates count data based on the overdispered Poisson distribution.
#' @param n The number of samples to generate.
#' @param network The underlying network (from create_network())
#' @param mu Mean expression for each gene. If length of mu is not equal to the number of
#' genes in the network, a random sample with replacement from mu will be used.
#' @param overdispersion A value > 0. Adjusts the amount of overdispersion in the generated counts.
#' @param intensity A value > 0. Used as the standard deviation of sampled edge weights.
#' @param k A value in [1, 2]. Adjusts the tail behavior of the distribution of generated counts.
#' @return A list containing the n by p matrix of samples, the underlying network,
#' and an n by p matrix containing mean expression values for each gene on each sample
#' after incorporating the underlying network.
#' @export
gen_gamma_poisson <- function(n, 
                              network, 
                              mu,
                              overdispersion = 1, 
                              intensity = 1,
                              k = 1.5) {
  if(k > 2 | k < 1) stop("k should be in [1, 2].")
  if(overdispersion <= 0) stop("overdispersion should be > 0")
  if(intensity <= 0) stop ("intensity should be > 0")
  
  p <- network$p
  if(length(mu) != p) {
    mu <- sample(mu, p, replace = TRUE)
  }
  max_mean_count <- max(mu)
  
  #S6: generate a vector of read counts for each sample.
  #x <- rpois(n * p, epsilon %*% t(mu * exp(beta %*% u)))
  x <- matrix(0, nrow = p, ncol = n) # Store samples in columns at first.
  mu <- matrix(mu, nrow = p, ncol = n + 1)
  adjust <- function(mu, edges) {
    x <- apply(edges, 2, function(edge) {
      index <- which(edge != 0)
      if(length(index) > 0) {
        return(mean(edge[index]))
      } else {
        return(0)
      }
    })
    val <- mu * exp(x)
    val <- ifelse(val <= 2 * max_mean_count, val, 2 * max_mean_count)
    return(val)
  }
  
  for(i in 1:n) {
    edges <- add_weight_to_network(network, intensity)$weight_matrix
    mu[, i + 1] <- adjust(mu[, 1], edges)
    theta <- rgamma(p,
                    shape = mu[, i + 1]^(2 - k) / overdispersion,
                    scale = mu[, i + 1]^(k - 1) * overdispersion)
    x[, i] <- rpois(p, theta)
  }
  
  x <- t(x) # Turn into n by p matrix.
  
  # Add column names to simulated dataset.
  if(is.null(network$node_names)) {
    colnames(x) <- paste(1:p)
  } else {
    colnames(x) <- network$node_names
  }
  
  if(any(is.na(x))) {
    warning(paste("NAs produced by `gen_gamma_poisson`: replacing",
                  sum(is.na(x)), "of", length(x), "values with 0"))
    x[which(is.na(x))] <- 0
  }
  
  return(list(x = x, 
              network = network, 
              mu = t(mu)))
}

attr(gen_gamma_poisson, "name") <- "gen_gamma_poisson"
