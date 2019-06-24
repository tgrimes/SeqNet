#' Get adjacency matrix
#' 
#' The adjacency matrix is constructed from all modules in a network.
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return An adjacency matrix with entry ij = 1 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @note The connections in an adjacency matrix and association matrix may differ
#' if the network contains multiple modules. The adjacency matrix only considers 
#' direct connections in the network, whereas the association matrix takes into 
#' account the fact that overlapping modules can create conditional dependencies
#' between two genes in seperate modules (i.e. genes that don't have a direct
#' connection in the graph).
#' @export
#' @examples 
#' # Create a random network with 10 nodes and add random edge weights.
#' nw <- random_network(10)
#' nw <- gen_partial_correlations(nw)
#' # Get adjacency matrix for the network or individual modules in the network.
#' get_adjacency_matrix(nw)
#' module <- nw$modules[[1]]
#' get_adjacency_matrix(module)
get_adjacency_matrix <- function(x, ...) {
  UseMethod("get_adjacency_matrix")
}

#' @inherit get_adjacency_matrix
#' @export
get_adjacency_matrix.default <- function(x, ...) {
  cat("get_adjacency_matrix() is defined for 'network' and 'network_module' objects.\n")
}

#' @inherit get_adjacency_matrix
#' @export
get_adjacency_matrix.matrix <- function(x, ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(dim(x)[1] != dim(x)[2])
    stop(paste0("'", deparse(substitute(x)), "' is not a square matrix."))
  if(!is_symmetric_cpp(x))   
    stop(paste0("'", deparse(substitute(x)), "' is not a symmetric matrix."))
  
  # Set any values that are near zero to exactly zero.
  x[abs(x) <= 10^-13] <- 0
  # Set all non-zero values to 1.
  x[x != 0] <- 1
  # Set diagonal values to zero.
  diag(x) <- 0
  return(x)
}



#' Get association matrix
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param tol A small tolerance threshold; any entry that is within `tol` from zero
#' is set to zero.
#' @param ... Additional arguments.
#' @return An association matrix with entry ij != 0 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @note The connections in an adjacency matrix and association matrix may differ
#' if the network contains multiple modules. The adjacency matrix only considers 
#' direct connections in the network, whereas the association matrix takes into 
#' account the fact that overlapping modules can create conditional dependencies
#' between two genes in seperate modules (i.e. genes that don't have a direct
#' connection in the graph).
#' @export
#' @examples 
#' # Create a random network with 10 nodes and add random edge weights.
#' nw <- random_network(10)
#' nw <- gen_partial_correlations(nw)
#' # Get adjacency matrix for the network or individual modules in the network.
#' get_association_matrix(nw)
#' module <- nw$modules[[1]]
#' get_association_matrix(module)
get_association_matrix <- function(x, tol = 10^-13, ...) {
  UseMethod("get_association_matrix")
}

#' @inherit get_association_matrix
#' @export
get_association_matrix.default <- function(x, tol = 10^-13, ...) {
  cat("get_association_matrix() is defined for 'network' and 'network_module' objects.\n")
}

#' @inherit get_association_matrix
#' @export
get_association_matrix.matrix <- function(x, tol = 10^-13, ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(dim(x)[1] != dim(x)[2])
    stop(paste0("'", deparse(substitute(x)), "' is not a square matrix."))
  if(!is_symmetric_cpp(x))   
    stop(paste0("'", deparse(substitute(x)),  "' is not a symmetric matrix."))
  if(!is_weighted(x)) 
    stop(paste0("'", deparse(substitute(x)), "' is not a weighted matrix."))

  diag(x) <- 0
  x[abs(x) < tol] <- 0 # Set entries near zero to zero.
  return(x)
}



#' Get the covariance matrix
#' 
#' The associations in each module are taken as partial correlations, and
#' the covariance matrix is calculated from these assuming that expression
#' for gene i is the weighted average over each module using 1/sqrt(m_i) 
#' as the weight, where m_i is the number of modules containing gene i.
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @note If a matrix is provided, it is assumed to be a partial correlation matrix;
#' a warning is given in this case. To avoid the warning message, convert the
#' matrix into a network object using 'create_network_from_association_matrix()'.
#' @return A covariance matrix.
#' @export
#' @examples 
#' # Create a random network with 10 nodes and add random edge weights.
#' nw <- random_network(10)
#' nw <- gen_partial_correlations(nw)
#' # Get covariance matrix for the network or individual modules in the network.
#' get_sigma(nw)
#' module <- nw$modules[[1]]
#' get_sigma(module)
get_sigma <- function(x, ...) {
  UseMethod("get_sigma")
}

#' @inherit get_sigma
#' @export
get_sigma.default <- function(x, ...) {
  # cat("get_sigma() is defined for 'network' and 'network_module' objects.\n")
}

#' @inherit get_sigma
#' @export
get_sigma.matrix <- function(x, ...) {
  if(any(diag(x) == 0)) {
    stop("Matrix 'x' should have nonzero entries along its diagonal.")
  }
  warning("Interpreting as a partial correlation matrix.")
  
  precision_matrix <- -x
  diag(precision_matrix) <- 1
  if(!is_PD(precision_matrix)) {
    stop(paste("The edge weights in the module do not correspond to a", 
               "positive definite precision matrix."))
  }
  sigma <- solve(precision_matrix)
  return(sigma)
}

#' Check if an object is weighted
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' object are weighted by 0s and 1s, and returns TRUE otherwise. If there
#' are no connections in the module, then this function returns TRUE.
#' @export
#' @examples 
#' # Create a random network with 10 nodes. 
#' nw <- random_network(10)
#' # The network, and hence all of its modules, are unweighted.
#' is_weighted(nw)
#' sapply(nw$modules, is_weighted)
#' # Add random weights to the connections.
#' nw <- gen_partial_correlations(nw)
#' # The network, and hence all of its modules, are now weighted.
#' is_weighted(nw)
#' sapply(nw$modules, is_weighted)
is_weighted <- function(x, ...) {
  UseMethod("is_weighted")
}

#' @inherit is_weighted
#' @export
is_weighted.default <- function(x, ...) {
  cat("is_weighted() is defined for 'network' and 'network_module' objects.\n")
}

#' @inherit is_weighted
#' @export
is_weighted.matrix <- function(x, ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(dim(x)[1] != dim(x)[2])
    stop(paste0("'", deparse(substitute(x)), "' is not a square matrix."))
  if(!is_symmetric_cpp(x))   
    stop(paste0("'", deparse(substitute(x)), "' is not a symmetric matrix."))
  
  x <- x[lower.tri(x)]
  
  if(any(!(x %in% c(0, 1)))) {
    return(TRUE)
  } else if(all(x == 0)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Removes the weights of all connections
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return The modified object.
#' @export
#' @examples 
#' # Create a random network with 10 nodes and add random edge weights.
#' nw <- random_network(10)
#' nw <- gen_partial_correlations(nw)
#' is_weighted(nw)
#' # Remove the edge weights from the network.
#' nw <- remove_weights(nw)
#' is_weighted(nw)
remove_weights <- function(x, ...) {
  UseMethod("remove_weights")
}

#' @inherit remove_weights
#' @export
remove_weights.default <- function(x, ...) {
  cat("remove_weights() is defined for 'network', 'network_module', and 'matrix' objects.\n")
}

#' @inherit remove_weights
#' @export
remove_weights.matrix <- function(x, ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(dim(x)[1] != dim(x)[2])
    stop(paste0("'", deparse(substitute(x)), "' is not a square matrix."))
  if(!is_symmetric_cpp(x))   
    stop(paste0("'", deparse(substitute(x)), "' is not a symmetric matrix."))
  
  x[x != 0] <- 1
  return(x)
}



#' Get node names
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return A vector containing the node names or node indicies.
#' @note Modules do not retain the names of each node, so the node indicies are
#' returned instead. These can be used to index into the vector of node
#' names obtained from the network.
#' @export
#' @examples 
#' # Create a random network with 10 nodes. 
#' nw <- random_network(10)
#' get_node_names(nw) # Default names are 1, 2, ..., 10.
#' nw <- set_node_names(nw, paste("node", 1:10, sep = "_"))
#' get_node_names(nw) # Print out updated node names.
#' # Modules only contain the indicies to nodes, not the node names
#' module <- nw$modules[[1]]
#' get_node_names(module)
#' # When converting the network to a matrix, node names appear as column names.
#' adj_matrix <- get_adjacency_matrix(nw)
#' colnames(adj_matrix) 
get_node_names <- function(x, ...) {
  UseMethod("get_node_names")
}

#' @inherit get_node_names
#' @export
get_node_names.default <- function(x, ...) {
  cat("get_node_names() is defined for 'network' and 'network_module' objects.\n")
}

#' @inherit get_node_names
#' @export
get_node_names.matrix <- function(x, ...) {
  node_names <- colnames(x)
  if(is.null(node_names)) {
    return(1:ncol(x))
  }
  return(node_names)
}




#' Rewire connections to a node
#' 
#' @param x The 'network', 'network_module', or 'matrix' object to modify.
#' @param node The node to rewire.
#' @param prob_rewire A value between 0 and 1, inclusive. Each connection to 'node' 
#' will be rewired with probability equal to 'prob_rewire'. Note, the degree of 
#' 'node' is unchanged after this operation.
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling nodes to rewire.
#' @param alpha A positive value used to parameterize the Beta distribution.
#' @param beta  A positive value used to parameterize the Beta distribution. 
#' @param epsilon A small constant added to the sampling probability of each node.
#' @param run_checks If TRUE and 'x' is a matrix, then it is checked that 'x' is an
#' adjacency matrix. This catches the case where 'x' is a weighted matrix, in which
#' case the weights are removed and a warning is given. 
#' @param ... Additional arguments.
#' @return The modified object.
#' @export
#' @examples 
#' # Create a random network with 10 nodes. 
#' nw <- random_network(10)
#' # Rewire connections to the first node.
#' nw_rewired <- rewire_connections_to_node(nw, 1)
#' # Plot the two networks for comparison
#' g <- plot(nw)
#' plot(nw_rewired, g) # Pass in g to mirror the layout.
#' # Or plot the differential network.
#' plot_network_diff(nw, nw_rewired, g)
rewire_connections_to_node <- function(x,
                                       node,
                                       prob_rewire = 1,
                                       weights = NULL,
                                       alpha = 100,
                                       beta = 1,
                                       epsilon = 10^-5,
                                       run_checks = TRUE,
                                       ...) {
  UseMethod("rewire_connections_to_node")
}

#' @inherit rewire_connections_to_node
#' @export
rewire_connections_to_node.default <- function(x,
                                               node,
                                               prob_rewire = 1,
                                               weights = NULL,
                                               alpha = 100,
                                               beta = 1,
                                               epsilon = 10^-5,
                                               run_checks = TRUE,
                                               ...) {
  cat("rewire_connections_to_node() is defined for 'network', 'network_module'",
      " and 'matrix' objects.\n")
}

#' @inherit rewire_connections_to_node
#' @export
rewire_connections_to_node.matrix <- function(x,
                                              node,
                                              prob_rewire = 1,
                                              weights = NULL,
                                              alpha = 100,
                                              beta = 1,
                                              epsilon = 10^-5,
                                              run_checks = TRUE,
                                              ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(run_checks && is_weighted(x)) {
    warning("'", deparse(substitute(x)), "' is weighted.",
            "Edge weights are removed prior to rewiring.")
    x <- remove_weights(x)
  }
  if(run_checks && !is_adjacency_cpp(x)) {
    # Stop if the matrix is not an adjacency matrix (even after removing weights).
    stop("'", deparse(substitute(x)), "' is not an adjacency matrix.")
  }
  
  nodes <- get_node_names(x)
  p <- length(nodes)
  
  if(node <= 0) 
    stop("Argument 'node' must be a positive integer.")
  if(prob_rewire < 0 || prob_rewire > 1) 
    stop("Argument 'prob_rewire' must be in [0, 1].")
  if(!is.null(weights) && (length(weights) != p || !is.numeric(weights))) 
    stop("'", deparse(substitute(weights)), "' is not a numeric vector of length", 
         " equal to the number of nodes in m.")
  
  if(is.null(weights)) {
    weights <- rep(0, p)
  }
  
  if(node %in% nodes) {
    node_index <- which(node == nodes)
    degree <- colSums(x)
    weights <- weights + degree
    
    # Rewire if the node has neighbors and is not connected to all other nodes.
    if(degree[node_index] > 0 && degree[node_index] != (p - 1)) {
      neighbors <- which(x[, node_index] != 0)
      available <- (1:p)[-c(node_index, neighbors)]
      
      # Determine how many connections to rewire
      n_rewire <- sum(runif(length(neighbors)) < prob_rewire)
      if(n_rewire > 0) {
        
        # Sample the connections to rewire.
        if(length(neighbors) == 1) {
          # Do not use sample() if only one neighbor.
          neighbors <- neighbors
        } else {
          neighbors <- sample(neighbors, n_rewire)
        }
        
        # Rewire each connection.
        for(old_neighbor in neighbors) {
          # Sample a new neighbor to wire to from the available nodes.
          if(length(available) == 1) {
            # Do not use sample() if only one available.
            new_neighbor <- available
          } else {
            # Sample high weight nodes with higher probability
            # new_neighbor <- sample(available, 1, prob = weights[available]^eta + epsilon)
            # new_neighbor <- sample(available, 1,
            #                        prob = ecdf_cpp(weights[available])^eta + epsilon)
            new_neighbor <- sample(available, 1,
                                   prob = pbeta(ecdf_cpp(weights[available]), alpha, beta) + epsilon)
            
          }
          # Rewire connection to new node.
          x[old_neighbor, node_index] <- 0
          x[node_index, old_neighbor] <- 0
          x[new_neighbor, node_index] <- 1
          x[node_index, new_neighbor] <- 1
          # The new neighbor is no longer available. Replace it with the old 
          # neighbor, which can now be selected.
          available[which(available == new_neighbor)] <- old_neighbor
          # Update weights.
          weights[new_neighbor] <- weights[new_neighbor] + 1
          weights[old_neighbor] <- weights[old_neighbor] - 1
        }
      }
    }
  } else {
    warning("'", deparse(substitute(node)), "' is not a node in 'x'.",
            "Returning 'x' unmodified.")
  }
  
  return(x)
}



#' Rewire connections
#' 
#' @param x The 'network', 'network_module', or 'matrix' object to modify.
#' @param prob_rewire A value between 0 and 1. The connections to each node 
#' will be rewired with probability equal to 'prob_rewire'. 
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling a node to rewire to.
#' @param alpha A positive value used to parameterize the Beta distribution.
#' @param beta  A positive value used to parameterize the Beta distribution. 
#' @param epsilon A small constant added to the sampling probability of each node.
#' @param run_checks If TRUE and 'x' is a matrix, then it is checked that 'x' is an
#' adjacency matrix. This catches the case where 'x' is a weighted matrix, in which
#' case the weights are removed and a warning is given. 
#' @param ... Additional arguments.
#' @return The modified module.
#' @note When applied to a network object, all modules in the network are
#' rewired. If
#' @export
#' @examples 
#' # Create a random network with 10 nodes. 
#' nw <- random_network(10)
#' # Rewire nodes in the network each with probability 1/2
#' nw_rewired <- rewire_connections(nw, 0.5)
#' # Plot the two networks for comparison
#' g <- plot(nw)
#' plot(nw_rewired, g) # Pass in g to mirror the layout.
#' # Or plot the differential network.
#' plot_network_diff(nw, nw_rewired, g)
rewire_connections <- function(x,
                               prob_rewire = 1,
                               weights = NULL,
                               alpha = 100,
                               beta = 1,
                               epsilon = 10^-5,
                               run_checks = TRUE,
                               ...) {
  UseMethod("rewire_connections")
}

#' @inherit rewire_connections
#' @export
rewire_connections.default <- function(x,
                                       prob_rewire = 1,
                                       weights = NULL,
                                       alpha = 100,
                                       beta = 1,
                                       epsilon = 10^-5,
                                       run_checks = TRUE,
                                       ...) {
  cat("rewire_connections() is defined for 'network_module'",
      " and 'matrix' objects.\n")
}

#' @inherit rewire_connections
#' @export
rewire_connections.matrix <- function(x,
                                      prob_rewire = 1,
                                      weights = NULL,
                                      alpha = 100,
                                      beta = 1,
                                      epsilon = 10^-5,
                                      run_checks = TRUE,
                                      ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(run_checks && is_weighted(x)) {
    warning("'", deparse(substitute(x)), "' is weighted.",
            "Edge weights are removed prior to rewiring.")
    x <- remove_weights(x)
  }
  if(run_checks && !is_adjacency_cpp(x)) {
    # Stop if the matrix is not an adjacency matrix (even after removing weights).
    stop("'", deparse(substitute(x)), "' is not an adjacency matrix.")
  }
  
  for(i in get_node_names(x)) {
    x <- rewire_connections_to_node(x, i, prob_rewire = prob_rewire,
                                    weights = weights, alpha = alpha, beta = beta, 
                                    epsilon = epsilon, ...)
  }
  
  return(x)
}



#' Remove connections to a node
#' 
#' @param x The 'network', 'network_module', or 'matrix' object to modify.
#' @param node The node to unwire.
#' @param prob_remove A value between 0 and 1. Each connection to 'node_index' 
#' will be removed with probability equal to 'prob_remove'.
#' @param run_checks If TRUE and 'x' is a matrix, then it is checked that 'x' is an
#' adjacency matrix. This catches the case where 'x' is a weighted matrix, in which
#' case the weights are removed and a warning is given. 
#' @param ... Additional arguments.
#' @return The modified adjacency matrix.
#' @export
#' @examples 
#' # Create a random network with 10 nodes. 
#' nw <- random_network(10)
#' # Remove all connections to node 1.
#' nw_rewired <- remove_connections_to_node(nw, 1, 1)
#' # Plot the two networks for comparison
#' g <- plot(nw)
#' plot(nw_rewired, g) # Pass in g to mirror the layout.
#' # Or plot the differential network.
#' plot_network_diff(nw, nw_rewired)
remove_connections_to_node <- function(x,
                                       node,
                                       prob_remove,
                                       run_checks = TRUE,
                                       ...) {
  UseMethod("remove_connections_to_node")
}

#' @inherit remove_connections_to_node
#' @export
remove_connections_to_node.default <- function(x,
                                               node,
                                               prob_remove,
                                               run_checks = TRUE,
                                               ...) {
  cat("remove_connections_to_node() is defined for 'network', 'network_module'",
      " and 'matrix' objects.\n")
}

#' @inherit remove_connections_to_node
#' @export
remove_connections_to_node.matrix <- function(x,
                                              node,
                                              prob_remove,
                                              run_checks = TRUE,
                                              ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(run_checks && is_weighted(x)) {
    warning("'", deparse(substitute(x)), "' is weighted.",
            "Edge weights are removed prior to rewiring.")
    x <- remove_weights(x)
  }
  if(run_checks && !is_adjacency_cpp(x)) {
    # Stop if the matrix is not an adjacency matrix (even after removing weights).
    stop("'", deparse(substitute(x)), "' is not an adjacency matrix.")
  }
  
  nodes <- get_node_names(x)
  p <- length(nodes)
  
  if(node <= 0) 
    stop("Argument 'node' must be a positive integer.")
  if(prob_remove < 0 || prob_remove > 1) 
    stop("Argument 'prob_remove' must be in [0, 1].")
  
  if(node %in% nodes) {
    node_index <- which(node == nodes)
    degree <- colSums(x)
    
    # Remove connections if the node has neighbors.
    if(degree[node_index] > 0) {
      neighbors <- which(x[, node_index] != 0)
      n_remove <- sum(runif(length(neighbors)) < prob_remove)
      
      if(n_remove > 0) {
        # Sample the connections to rewire.
        if(length(neighbors) == 1) {
          # Do not use sample() if only one neighbor.
          neighbors <- neighbors
        } else {
          neighbors <- sample(neighbors, n_remove)
        }
        
        # Remove connections
        x[neighbors, node_index] <- 0
        x[node_index, neighbors] <- 0
      }
    }
  } else {
    warning("'", deparse(substitute(node)), "' is not a node in 'x'.",
            "Returning 'x' unmodified.")
  }
  
  return(x)
}


#' Remove connections in a network
#' 
#' @param x The 'network', 'network_module', or 'matrix' object to modify.
#' @param prob_remove A value between 0 and 1. Each edge will be removed with 
#' probability equal to 'prob_remove'.
#' @param run_checks If TRUE and 'x' is a matrix, then it is checked that 'x' is an
#' adjacency matrix. This catches the case where 'x' is a weighted matrix, in which
#' case the weights are removed and a warning is given. 
#' @param ... Additional arguments.
#' @return The modified adjacency matrix.
#' @export
#' @examples 
#' # Create a random network with 10 nodes. 
#' nw <- random_network(20)
#' # Remove connections in the network each with probability 1/2.
#' nw_rewired <- remove_connections(nw, 0.5)
#' # Plot the two networks for comparison
#' g <- plot(nw)
#' plot(nw_rewired, g) # Pass in g to mirror the layout.
#' # Or plot the differential network.
#' plot_network_diff(nw, nw_rewired)
remove_connections <- function(x,
                               prob_remove,
                               run_checks = TRUE,
                               ...) {
  UseMethod("remove_connections")
}

#' @inherit remove_connections
#' @export
remove_connections.default <- function(x,
                                       prob_remove,
                                       run_checks = TRUE,
                                       ...) {
  cat("remove_connections() is defined for ",
      "'matrix' objects.\n")
}

#' @inherit remove_connections
#' @export
remove_connections.matrix <- function(x,
                                      prob_remove,
                                      run_checks = TRUE,
                                      ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(run_checks && is_weighted(x)) {
    warning("'", deparse(substitute(x)), "' is weighted.",
            "Edge weights are removed prior to rewiring.")
    x <- remove_weights(x)
  }
  if(run_checks && !is_adjacency_cpp(x)) {
    # Stop if the matrix is not an adjacency matrix (even after removing weights).
    stop("'", deparse(substitute(x)), "' is not an adjacency matrix.")
  }
  
  p <- ncol(x)
  
  if(prob_remove < 0 || prob_remove > 1) 
    stop("Argument 'prob_remove' must be in [0, 1].")
  
  # Obtain list of edges in the network.
  edges <- edges_from_adjacency_cpp(x)
  
  n_remove <- sum(runif(nrow(edges)) < prob_remove)
  if(n_remove > 0) {
    edge_index <- sample(1:nrow(edges), n_remove)
    
    # Remove connections
    if(length(edge_index) == 1) {
      x[edges[edge_index, 1], edges[edge_index, 2]] <- 0
      x[edges[edge_index, 2], edges[edge_index, 1]] <- 0
    } else {
      x[edges[edge_index, 1:2]] <- 0
      x[edges[edge_index, 2:1]] <- 0
    }

  }
  
  return(x)
}
