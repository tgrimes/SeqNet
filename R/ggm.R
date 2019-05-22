#' Generate observations from a Gaussian graphical model.
#' 
#' Generates data based on the multivariate normal distribution parameterized by
#' a zero mean vector and a covariance matrix. Observations are generated for
#' each module in the network individually, and the covariance matrix is set to
#' the inverse of the standardized association matrix for the module. 
#' Observations are combined for gene i by taking the sum across the m_i modules 
#' containing it and dividing by sqrt(m_i). 
#' @param n The number of samples to generate. If multiple networks are provided,
#' n samples are generated per network.
#' @param ... The 'network' object(s) to generate data from. Can be a single
#' network, many networks, or a single list of networks.
#' @return A list containing the n by p matrix of samples and the 'network'
#' object used to generate them.
#' @export
gen_gaussian <- function(n, ...) {
  if(n <= 0) {
    stop("n must be positive.")
  }
  
  # Handle network arguments
  network_list <- get_network_arguments(...)
  
  # If a single, unlisted network was provided, then the generated results are 
  # returned in the same format, unlisted.
  if(length(list(...)) == 1 && class(list(...)[[1]]) == "network") {
    single_network <- TRUE
  } else {
    single_network <- FALSE
  }
  
  if(!all(sapply(network_list, is_weighted))) {
    warning(paste("Network argument(s) are not all weighted.",
                  "Using gen_partial_correlations() to create weights."))
    weighted_networks <- gen_partial_correlations(network_list)
    network_list <- weighted_networks
  }
  
  n_networks <- length(network_list)
  x_list <- vector("list", n_networks)
  m_list <- vector("list", n_networks)
  for(i in 1:n_networks) {
    # Generated observations for the network.
    x <- matrix(0, n, network_list[[i]]$p)
    # Number of modules containing each node - updated after each iteration below.
    m <- rep(0, network_list[[i]]$p)
    
    if(length(network_list[[i]]$modules) > 0) {
      for(j in 1:length(network_list[[i]]$modules)) {
        nodes <- network_list[[i]]$modules[[j]]$nodes
        sigma <- get_sigma(network_list[[i]]$modules[[j]])
        
        x[, nodes] <- x[, nodes] + mvtnorm::rmvnorm(n, sigma = sigma)
        m[nodes] <- m[nodes] + 1
      }
    } else {
      x <- mvtnorm::rmvnorm(n, sigma = diag(1, ncol(x)))
      m <- m + 1
    }
    
    
    # Generate observations for nodes not in any module
    index <- which(m == 0)
    x[, index] <- rnorm(n * length(index))
    m[index] <- 1
    
    # Divide column i by 1 / sqrt(m_i) for i = 1, ..., p.
    x <- x %*% diag(1 / sqrt(m))
    
    # Standardize columns by dividing by standard deviation.
    x <- x %*% diag(1 / sqrt(diag(get_sigma(network_list[[i]]))))
    
    if(any(is.na(x))) {
      stop(paste("NAs produced in gen_gaussian() for network", i, "module", j))
    }
    
    x_list[[i]] <- x
    
    # Add column names to simulated dataset.
    colnames(x_list[[i]]) <- network_list[[i]]$node_names
  }
  
  result <- list(x = x_list,
                 network = network_list)
  
  # If a single network was provided, unlist the elements in the result.
  if(single_network) {
    result$x <- result$x[[1]]
    result$network <- result$network[[1]]
  } else {
    # Use network names to label each generated dataset.
    names(result$x) <- names(network_list)
  }
  
  return(result)
}








#' Generate partial correlations for a list of networks.
#' 
#' Random partial correlations are generated to weigh the network connections. 
#' If multiple networks are provided, the networks must contain the same nodes
#' and the same modules (the connections within modules may differ). Any 
#' connection that is common across different networks will also have the same 
#' partial correlation weight across networks.
#' @param ... The 'network' objects to modify.
#' @param k An integer that ensures the matrix inverse is numerically stable. 
#' k = 2.5 is default; higher values will allow for larger values of
#' partial correlations (and will result in a wider distribution of 
#' Pearson correlations).
#' @param rweights A generator for initial weights in the network. By default, 
#' values are generated uniformly from (-1, -0.5) U (0.5, 1). The weights will
#' be adjusted so that the sign of a generated weight and the sign of the
#' corresponding partial correlation agree.
#' @return An updated network object containing random weights. If multiple
#' networks were provided, then a list of network objects is returned.
#' @export
gen_partial_correlations <- function(...,
                                     k = 2.5,
                                     rweights = function(n) (-1)^rbinom(n, 1, 0.5) * runif(n, 0.5, 1)) {
  # Check 'k'.
  if(k < 1) 
    stop("Argument 'k' must be >= 1.")
  
  # Handle network arguments
  network_list <- get_network_arguments(...)
  
  # If a single network was provided, the updated network is returned unlisted.
  if(length(list(...)) == 1 && class(list(...)[[1]]) == "network") {
    single_network <- TRUE
  } else {
    single_network <- FALSE
  }
  
  # At this point a list of network objects is available.
  # Next steps are to
  #   0. Check that the networks contain the same nodes,
  #   1. initialize random association weights for each possible connection,
  #   2. find an adjustment to ensure invertibility of association matrix, and
  #   3. use the maximum adjustment to obtain precision matricies.
  
  if(!all_networks_contain_same_nodes(network_list)) {
    stop(paste("The networks do not contain the same nodes.",
               "The node names and order must be identical for each network."))
  }
  if(!all_networks_contain_same_modules(network_list)) {
    stop(paste("The networks do not contain the same modules."))
  }
  
  n_networks <- length(network_list)
  n_modules <- length(network_list[[1]]$modules)
  
  if(n_modules > 0) {
    for(i in 1:n_modules) {
      module_list <- lapply(network_list, function(network) network$modules[[i]])
      node_names <- module_list[[1]]$nodes
      p <- length(node_names)
      
      # Generate association weights for every possible connection in the module.
      # Associations along diagonal are zero; these are adjusted later.
      weights <- matrix(0, p, p)
      m <- p * (p - 1) / 2
      # Use negative weight so that partial correlations have same sign as
      # the generated weight.
      weights[lower.tri(weights)] <- -rweights(m)
      weights <- weights + t(weights)
      
      # Obtain an association matrix for each network using the generated values.
      weight_matrix_list <- lapply(module_list, function(module) {
        weight_matrix <- weights * get_adjacency_matrix(module)
        weight_matrix
      })
      
      # Ensure each association matrix is invertible by adjusting the diagonal.
      eigen_val_list <- lapply(weight_matrix_list, function(m) eigen(m)$values)
      adjustment_list <- lapply(eigen_val_list, function(lambda) {
        (max(lambda) * 10^-k - min(lambda))
      })
      
      # Find the maximum adjustment and apply to each association matrix.
      # Each matrix is now invertible and is interpreted as the precision matrix
      adjustment <- diag(max(unlist(adjustment_list)), p)
      precision_list <- lapply(weight_matrix_list, function(m) m + adjustment)
      
      # Standardize with 1s along diagonal.
      # Note: the precision matrix is related to negative partial correlations.
      # pcor_ij = -precision_ij/sqrt(precision_ii*precision_jj), for i != j.
      pcor_matrix_list <- lapply(precision_list, function(Omega) {
        if(all(Omega[lower.tri(Omega)] == 0)) {
          # If the module has no connections, set diagonal matrix.
          pcor_matrix <- diag(1, nrow(Omega))
        } else {
          pcor_matrix <- -cov2cor(Omega)
          diag(pcor_matrix) <- 1 # -diag(-cov2cor(Omega)) = 1
        }
        colnames(pcor_matrix) <- node_names
        return(pcor_matrix)
      })
      
      for(j in 1:n_networks) {
        # If the module contains edges, then set the edge weights.
        if(!is.null(network_list[[j]]$modules[[i]]$edges)) {
          network_list[[j]]$modules[[i]] <- 
            set_module_weights(network_list[[j]]$modules[[i]],
                               pcor_matrix_list[[j]])
        }
      }
    }
  }
  
  # If a single network was provided, return the updated network unlisted.
  if(single_network) {
    network_list <- network_list[[1]]
  }
  return(network_list)
}


#' Internal function to check if a list of networks all contain the same nodes.
#' 
#' @param network_list A list of 'network' objects.
#' @return A logical value; TRUE indicates the networks contain the same nodes,
#' FALSE indicates otherwise.
all_networks_contain_same_nodes <- function(network_list) {
  n_networks <- length(network_list)
  if(n_networks > 1) {
    nodes <- network_list[[1]]$node_names
    for(i in 2:n_networks) {
      # Check that networks have equal number of nodes and identical node names.
      temp <- network_list[[i]]$node_names
      if(length(nodes) != length(temp)) {
        return(FALSE)
      } 
      if(!all(nodes == temp)) {
        return(FALSE)
      }
    }
  }
  
  return(TRUE)
}


#' Internal function used to extract 'network' objects from argument list.
#' 
#' @param ... The 'network' object(s) or list of networks.
#' @return A list of 'network' objects.
get_network_arguments <- function(...) {
  network_list <- list(...)
  if(length(network_list) == 0) {
    stop("At least one 'network' object must be provided.")
  }
  
  # If a list of networks were provided, unlist to obtain original list.
  if(class(network_list[[1]]) == "list") {
    network_list <- network_list[[1]]
    network_names <- colnames(network_list)
  } else {
    # Otherwise, keep list and get the network variable names
    network_names <- sapply(substitute(list(...)), deparse)[-1]
    if(!is.null(names(network_list))) {
      # If any arguments are named, use those names.
      temp <- names(network_list)
      index <- which(!is.na(temp) & temp != "")
      network_names[index] <- temp[index]
    }
    names(network_list) <- network_names
  }
  
  # Check that network_list contains all 'network' objects.
  index <- which(sapply(network_list, class) != "network")
  if(length(index) == 1) {
    stop(paste0("Argument '", 
                network_names[index], 
                "' is not a 'network' object."))
  } else if(length(index) > 1) {
    stop(paste0("Arguments ",
                paste(paste0("'", network_names[index], "'"), collapse = ", "), 
                " are not 'network' objects."))
  }
  
  return(network_list)
}

#' Internal function to check if a list of networks all contain the same modules.
#' 
#' @param network_list A list of 'network' objects.
#' @return A logical value; TRUE indicates the networks contain the same modules,
#' FALSE indicates otherwise. Note, this only checks that the modules contain
#' the same nodes - the structure of the modules are allowed to differ.
all_networks_contain_same_modules <- function(network_list) {
  n_networks <- length(network_list)
  if(n_networks > 1) {
    n_modules <- length(network_list[[1]]$modules)
    if(n_modules > 0) {
      modules <- network_list[[1]]$modules
      for(i in 1:n_modules) {
        for(j in 2:n_networks) {
          temp <- network_list[[j]]$modules[[i]]
          if(length(modules[[i]]$nodes) != length(temp$nodes)) {
            return(FALSE)
          } 
          if(!all(modules[[i]]$nodes == temp$nodes)) {
            return(FALSE)
          }
        }
      }
    }
  }
  
  return(TRUE)
}