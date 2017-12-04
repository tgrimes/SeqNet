# Gamma-Poisson model

# Requires: 
#   SimSeq - for default reference dataset
library(SimSeq)

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
#' @param overdis_param The overdispersion parameters (shape, scale) used 
#'   in the Gamma distribution.
#' @param intensity_param The standard deviation of sampled edge weights.
#' @param k Adjusts the tail behavior in the simulated expression distribution.
#' @return A list containing the n by p matrix of samples and the true
#'   association matrix.
#' @export
gen_gamma_poisson <- function(n, network, mu,
                              overdispersion_param = 1, 
                              intensity_param = 1,
                              k = 1) {
  p <- network$p
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
    val <- ifelse(val <= max_mean_count, val, 2 * mu)
    return(val)
  }
  
  for(i in 1:n) {
    edges <- add_weight_to_network(network, intensity_param)$weight_matrix
    mu[, i + 1] <- adjust(mu[, 1], edges)
    theta <- rgamma(p,
                    shape = mu[, i + 1]^(2 - k) / overdispersion_param,
                    scale = mu[, i + 1]^(k - 1) * overdispersion_param)
    x[, i] <- rpois(p, theta)
  }
  
  x <- t(x) # Turn into n by p matrix.
  colnames(x) <- colnames(network$adj_matrix)
  
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
