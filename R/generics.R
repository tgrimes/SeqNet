#' Get adjacency matrix
#' 
#' The adjacency matrix is constructed from all modules in a network.
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return An adjacency matrix with entry ij = 1 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_adjacency_matrix <- function(x, ...) {
  UseMethod("get_adjacency_matrix")
}

#' Get adjacency matrix
#' 
#' The adjacency matrix is constructed from all modules in a network.
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return An adjacency matrix with entry ij = 1 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_adjacency_matrix.default <- function(x, ...) {
  cat("get_adjacency_matrix() is defined for 'network' and 'network_module' objects.\n")
}

#' Get adjacency matrix
#' 
#' All off-diagonal entries that are near zero are set to zero. All remaining 
#' non-zero values are set to 1, and the diagonal is set to zero.
#' @param x A 'matrix' object.
#' @param ... Additional arguments.
#' @return An adjacency matrix with entry ij = 1 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_adjacency_matrix.matrix <- function(x, ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(dim(x)[1] != dim(x)[2])
    stop(paste0("'", deparse(substitute(x)), "' is not a square matrix."))
  if(max(abs(x - t(x))) > 10^-13)   
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
#' @param ... Additional arguments.
#' @return An association matrix with entry ij != 0 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_association_matrix <- function(x, ...) {
  UseMethod("get_association_matrix")
}

#' Get association matrix
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' #' @return An association matrix with entry ij != 0 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_association_matrix.default <- function(x, ...) {
  cat("get_association_matrix() is defined for 'network' and 'network_module' objects.\n")
}

#' Get association matrix
#' 
#' All off-diagonal values in the matrix are left unchanged. The diagonal
#' is set to zero.
#' @param x A 'matrix' object.
#' @param tol A numeric value. Associations within tol from zero are set to zero.
#' @param ... Additional arguments.
#' @return Returns the matrix with diagonal elements set to zero.
#' @export
get_association_matrix.matrix <- function(x, tol = 10^-13, ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(dim(x)[1] != dim(x)[2])
    stop(paste0("'", deparse(substitute(x)), "' is not a square matrix."))
  if(max(abs(x - t(x))) > 10^-13)   
    stop(paste0("'", deparse(substitute(x)),  "' is not a symmetric matrix."))
  if(!is_weighted(x)) 
    stop(paste0("'", deparse(substitute(x)), "' is not a weighted matrix."))
  
  diag(x) <- 0
  x[abs(x) < tol] <- 0
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
#' @note If a matrix is provided, it is assumed to be a partial correlation matrix.
#' A warning is given in this case. To avoid the warning message, convert the
#' matrix into a network object using 'create_network_from_association_matrix()'.
#' @return A covariance matrix.
#' @export
get_sigma <- function(x, ...) {
  UseMethod("get_sigma")
}

#' Get the covariance matrix
#' 
#' The associations in each module are taken as partial correlations, and
#' the covariance matrix is calculated from these assuming that expression
#' for gene i is the weighted average over each module using 1/sqrt(m_i) 
#' as the weight, where m_i is the number of modules containing gene i.
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @note If a matrix is provided, it is assumed to be a partial correlation matrix.
#' A warning is given in this case. To avoid the warning message, convert the
#' matrix into a network object using 'create_network_from_association_matrix()'.
#' @return A covariance matrix.
#' @export
get_sigma.default <- function(x, ...) {
  cat("get_sigma() is defined for 'network' and 'network_module' objects.\n")
}

#' Get the covariance matrix
#' 
#' The associations in each module are taken as partial correlations, and
#' the covariance matrix is calculated from these assuming that expression
#' for gene i is the weighted average over each module using 1/sqrt(m_i) 
#' as the weight, where m_i is the number of modules containing gene i.
#' @param x A 'matrix' object.
#' @param ... Additional arguments.
#' @note The matrix is assumed to be a partial correlation matrix.
#' A warning is given in this case. To avoid the warning message, convert the
#' matrix into a network object using 'create_network_from_association_matrix()'.
#' @return A covariance matrix.
#' @export
get_sigma.matrix <- function(x, ...) {
  if(any(diag(x) == 0)) {
    stop("Matrix 'x' should have nonzero entries along its diagonal.")
  }
  warning("Interpreting as a partial correlation matrix.")
  
  precision_matrix <- -x
  if(all(precision_matrix == 0)) {
    # If there are no connections in x, return the identify matrix.
    return(diag(1, nrow(precision_matrix)))
  }
  diag(precision_matrix) <- 1
  
  if(any(diag(chol(precision_matrix)) < 0)) {
    warning(paste("The edge weights in the module do not correspond to a", 
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
is_weighted <- function(x, ...) {
  UseMethod("is_weighted")
}

#' Check if an object is weighted
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return Returns FALSE if all of the connections in the 
#' object are weighted by 0s and 1s, and returns TRUE otherwise. If there
#' are no connections in the module, then this function returns TRUE.
#' @export
is_weighted.default <- function(x, ...) {
  cat("is_weighted() is defined for 'network' and 'network_module' objects.\n")
}

#' Check if an object is weighted
#' 
#' @param x A 'matrix' object.
#' @param ... Additional arguments.
#' @return Returns FALSE if all of entries in the matrix are either zero or one, 
#' and returns TRUE otherwise. However, if there are no connections, i.e. all
#' off-diagonal entries are zero, then this function returns TRUE.
#' @export
is_weighted.matrix <- function(x, ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(dim(x)[1] != dim(x)[2])
    stop(paste0("'", deparse(substitute(x)), "' is not a square matrix."))
  if(max(abs(x - t(x))) > 10^-13)   
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
remove_weights <- function(x, ...) {
  UseMethod("remove_weights")
}

#' Removes the weights of all connections
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return The modified object.
#' @export
remove_weights.default <- function(x, ...) {
  cat("remove_weights() is defined for 'network', 'network_module', and 'matrix' objects.\n")
}

#' Removes the weights of all connections
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return The modified object.
#' @export
remove_weights.matrix <- function(x, ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  if(dim(x)[1] != dim(x)[2])
    stop(paste0("'", deparse(substitute(x)), "' is not a square matrix."))
  if(max(abs(x - t(x))) > 10^-13)   
    stop(paste0("'", deparse(substitute(x)), "' is not a symmetric matrix."))
  
  x[x != 0] <- 1
  return(x)
}



#' Get node names
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return A vector containing the node names or node indicies.
#' @export
get_node_names <- function(x, ...) {
  UseMethod("get_node_names")
}

#' Get node names
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return A vector containing the node names or node indicies.
#' @export
get_node_names.default <- function(x, ...) {
  cat("get_node_names() is defined for 'network' and 'network_module' objects.\n")
}

#' Get node names
#' 
#' @param x A 'matrix' object.
#' @param ... Additional arguments.
#' @return A vector containing the node names or node indicies.
#' @export
get_node_names.matrix <- function(x, ...) {
  node_names <- colnames(x)
  if(is.null(node_names)) {
    return(1:ncol(x))
  }
  return(node_names)
}




#' Rewire connections to a node.
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param node The node to rewire.
#' @param ... Additional arguments.
#' @return The modified object.
#' @export
rewire_connections_to_node <- function(x, node, ...) {
  UseMethod("rewire_connections_to_node")
}

#' Rewire connections to a node.
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param node The node to rewire.
#' @param ... Additional arguments. 
#' @return The modified object.
#' @export
rewire_connections_to_node.default <- function(x, node, ...) {
  cat("rewire_connections_to_node() is defined for 'network', 'network_module'",
      " and 'matrix' objects.\n")
}

#' Rewire connections to a node.
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param node The node to rewire.
#' @param prob_rewire A value between 0 and 1. Each connection to 'node' 
#' will be rewired with probability equal to 'prob_rewire'. Note, the degree of 
#' 'node' is unchanged after this operation.
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling nodes to rewire.
#' @param exponent The exponent used for weighted sampling. When exponent = 0,
#' nodes are sampled uniformly. When exponent > 0, the sampling probability
#' is based on node weights.
#' @param ... Additional arguments.
#' @return The modified object x.
#' @export
rewire_connections_to_node.matrix <- function(x,
                                              node,
                                              prob_rewire,
                                              weights = NULL,
                                              exponent = 0,
                                              ...) {
  if(!(class(x) == "matrix")) 
    stop(paste0("'", deparse(substitute(x)), "' is not a 'matrix' object."))
  
  nodes <- get_node_names(x)
  p <- length(nodes)
  
  if(node <= 0) 
    stop("Argument 'node' must be a positive integer.")
  if(prob_rewire < 0 || prob_rewire > 1) 
    stop("Argument 'prob_rewire' must be in [0, 1].")
  if(!is.null(weights) && (length(weights) != p || !is.numeric(weights))) 
    stop("'", deparse(substitute(weights)), "' is not a numeric vector of length", 
         " equal to the number of nodes in m.")
  if(length(exponent) != 1 || !is.numeric(exponent)) 
    stop("'", deparse(substitute(exponent)), "' is not positive numeric value.")
  
  
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
          # Sample low weight neighbors with higher probability.
          neighbors <- sample(neighbors, n_rewire,
                              prob = (1 - ecdf_cpp(weights[neighbors]))^exponent + 10^-16)
        }
        
        # Rewire each connection.
        for(old_neighbor in neighbors) {
          # Sample a new neighbor to wire to from the available nodes.
          if(length(available) == 1) {
            # Do not use sample() if only one available.
            new_neighbor <- available
          } else {
            # Sample high weight nodes with higher probability
            new_neighbor <- sample(available, 1,
                                   prob = ecdf_cpp(weights[available])^exponent + 10^-16)
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



#' Rewire connections.
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments.
#' @return The modified object.
#' @export
rewire_connections <- function(x, ...) {
  UseMethod("rewire_connections")
}

#' Rewire connections.
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param ... Additional arguments. 
#' @return The modified object.
#' @export
rewire_connections.default <- function(x, ...) {
  cat("rewire_connections() is defined for 'network', 'network_module'",
      " and 'matrix' objects.\n")
}

#' Rewire connections.
#' 
#' @param x Either a 'network', 'network_module', or 'matrix' object.
#' @param prob_rewire A value between 0 and 1. The connections to each node 
#' will be rewired with probability equal to 'prob_rewire'. 
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling a node to rewire to.
#' @param exponent The exponent used for weighted sampling. When exponent = 0,
#' nodes are sampled uniformly. When exponent > 0, the sampling probability
#' is based on node weights.
#' @param ... Additional arguments.
#' @return The modified object x.
#' @export
rewire_connections.matrix <- function(x,
                                      prob_rewire,
                                      weights = NULL,
                                      exponent = 0,
                                      ...) {
  for(i in 1:ncol(x)) {
    x <- rewire_connections_to_node(x, i, prob_rewire = prob_rewire,
                                    weights = weights, exponent = exponent, ...)
  }
  
  return(x)
}



#' Remove connections to a node.
#' 
#' @param x A 'network', 'network_module', or 'matrix' object to modify. 
#' @param node The node to rewire.
#' @param ... Additional parameters.
#' @return The modified object.
#' @export
remove_connections_to_node <- function(x, node, ...) {
  UseMethod("remove_connections_to_node")
}

#' Remove connections to a node.
#' 
#' @param x A 'network', 'network_module', or 'matrix' object to modify. 
#' @param node The node to rewire.
#' @param ... Additional parameters.
#' @return The modified object.
#' @export
remove_connections_to_node.default <- function(x, node, ...) {
  cat("remove_connections_to_node() is defined for 'network', 'network_module'",
      " and 'matrix' objects.\n")
}

#' Remove connections to a node.
#' 
#' @param x An adjacency matrix to modify.
#' @param node The node to unwire.
#' @param prob_remove A value between 0 and 1. Each connection to 'node_index' 
#' will be removed with probability equal to 'prob_remove'.
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling neighbors to unwire from.
#' @param exponent The exponent used for weighted sampling. When exponent = 0,
#' neighboring nodes are sampled uniformly. When exponent > 0, the sampling 
#' probability is based on node weights.
#' @param ... Additional arguments.
#' @return The modified adjacency matrix.
#' @export
remove_connections_to_node.matrix <- function(x,
                                              node,
                                              prob_remove,
                                              weights = NULL,
                                              exponent = 0,
                                              ...) {
  if(!(class(x) == "matrix")) 
    stop("'", deparse(substitute(x)), "' is not a 'matrix' object.")
  
  nodes <- get_node_names(x)
  p <- length(nodes)
  
  if(node <= 0) 
    stop("Argument 'node' must be a positive integer.")
  if(prob_remove < 0 || prob_remove > 1) 
    stop("Argument 'prob_remove' must be in [0, 1].")
  if(!is.null(weights) && (length(weights) != p || !is.numeric(weights))) 
    stop("'", deparse(substitute(weights)), "' is not a numeric vector of length", 
         " equal to the number of nodes in 'x'.")
  if(length(exponent) != 1 || !is.numeric(exponent)) 
    stop("'", deparse(substitute(exponent)), "' is not positive numeric value.")
  
  
  
  if(is.null(weights)) {
    weights <- rep(0, p)
  }
  
  if(node %in% nodes) {
    node_index <- which(node == nodes)
    degree <- colSums(x)
    weights <- weights + degree
    
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
          # Sample low weight neighbors with higher probability.
          neighbors <- sample(neighbors, n_remove,
                              prob = (1 - ecdf_cpp(weights[neighbors]))^exponent + 10^-16)
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


#' Remove connections in a network.
#' 
#' @param x A 'matrix' object to modify. 
#' @param ... Additional parameters.
#' @return The modified object.
#' @export
remove_connections <- function(x, ...) {
  UseMethod("remove_connections")
}

#' Remove connections in a network.
#' 
#' @param x A 'matrix' object to modify. 
#' @param ... Additional parameters.
#' @return The modified object.
#' @export
remove_connections.default <- function(x, ...) {
  cat("remove_connections() is defined for ",
      "'matrix' objects.\n")
}

#' Remove connections in a network.
#' 
#' @param x An adjacency matrix to modify.
#' @param prob_remove A value between 0 and 1. Each edge will be removed with 
#' probability equal to 'prob_remove'.
#' @param weights (Optional) A vector of weights for each node. These are used,
#' in addition to the degree of each node, when sampling edges to remove.
#' @param exponent The exponent used for weighted sampling. When exponent = 0,
#' edges are sampled uniformly. When exponent > 0, the sampling probability is
#' based on the connected node weights.
#' @param ... Additional arguments.
#' @return The modified adjacency matrix.
#' @export
remove_connections.matrix <- function(x,
                                      prob_remove,
                                      weights = NULL,
                                      exponent = 0,
                                      ...) {
  if(!(class(x) == "matrix")) 
    stop("'", deparse(substitute(x)), "' is not a 'matrix' object.")
  
  p <- ncol(x)
  
  if(prob_remove < 0 || prob_remove > 1) 
    stop("Argument 'prob_remove' must be in [0, 1].")
  if(!is.null(weights) && (length(weights) != p || !is.numeric(weights))) 
    stop("'", deparse(substitute(weights)), "' is not a numeric vector of length", 
         " equal to the number of nodes in m.")
  if(length(exponent) != 1 || !is.numeric(exponent)) 
    stop("'", deparse(substitute(exponent)), "' is not positive numeric value.")
  
  
  if(is.null(weights)) {
    weights <- rep(0, p)
  }
  weights <- weights + colSums(x)
  weights <- unname(weights)
  
  # Obtain list of edges in the network.
  edges <- edges_from_adjacency_cpp(x)
  weights <- weights[edges[, 1]] + weights[edges[, 2]]
  
  n_remove <- sum(runif(nrow(edges)) < prob_remove)
  
  if(n_remove > 0) {
    
    # Sample high weight edges with higher probability.
    edge_index <- sample(1:nrow(edges), n_remove,
                         prob = (ecdf_cpp(weights))^exponent + 10^-16)
    
    # Remove connections
    x[matrix(edges[edge_index, ], ncol = 2)] <- 0
    x[matrix(edges[edge_index, 2:1], ncol = 2)] <- 0
  }
  
  return(x)
}
