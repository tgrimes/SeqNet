#' Get adjacency matrix
#' 
#' The adjacency matrix is constructed from all modules in a network.
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return An adjacency matrix with entry ij = 1 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_adjacency_matrix <- function(...) {
  UseMethod("get_adjacency_matrix")
}

#' Get adjacency matrix
#' 
#' The adjacency matrix is constructed from all modules in a network.
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return An adjacency matrix with entry ij = 1 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_adjacency_matrix.default <- function(...) {
  cat("get_adjacency_matrix() is defined for 'network' and 'network_module' objects.\n")
}

#' Get adjacency matrix
#' 
#' All off-diagonal entries that are near zero are set to zero. All remaining 
#' non-zero values are set to 1, and the diagonal is set to zero.
#' @param m A 'matrix' object.
#' @return An adjacency matrix with entry ij = 1 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_adjacency_matrix.matrix <- function(m, ...) {
  if(!(class(m) == "matrix")) 
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a 'matrix' object."))
  if(dim(m)[1] != dim(m)[2])
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a square matrix."))
  if(max(abs(m - t(m))) > 10^-13)   
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a symmetric matrix."))
  
  # Set any values that are near zero to exactly zero.
  m[abs(m) <= 10^-13] <- 0
  # Set all non-zero values to 1.
  m[m != 0] <- 1
  # Set diagonal values to zero.
  diag(m) <- 0
  return(m)
}



#' Get association matrix
#' 
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return An association matrix with entry ij != 0 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_association_matrix <- function(...) {
  UseMethod("get_association_matrix")
}

#' Get association matrix
#' 
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return An association matrix with entry ij != 0 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_association_matrix.default <- function(...) {
  cat("get_association_matrix() is defined for 'network' and 'network_module' objects.\n")
}

#' Get association matrix
#' 
#' All off-diagonal values in the matrix are left unchanged. The diagonal
#' is set to zero.
#' @param m A 'matrix' object.
#' @return Returns the matrix with diagonal elements set to zero.
#' @export
get_association_matrix.matrix <- function(m, ...) {
  if(!(class(m) == "matrix")) 
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a 'matrix' object."))
  if(dim(m)[1] != dim(m)[2])
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a square matrix."))
  if(max(abs(m - t(m))) > 10^-13)   
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a symmetric matrix."))
  if(!is_weighted(m)) 
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a weighted matrix."))
  
  diag(m) <- 0
  return(m)
}



#' Get the covariance matrix
#' 
#' The associations in each module are taken as partial correlations, and
#' the covariance matrix is calculated from these assuming that expression
#' for gene i is the weighted average over each module using 1/sqrt(m_i) 
#' as the weight, where m_i is the number of modules containing gene i.
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @note If a matrix is provided, it is assumed to be a partial correlation matrix.
#' A warning is given in this case. To avoid the warning message, convert the
#' matrix into a network object using 'create_network_from_association_matrix()'.
#' @return A covariance matrix.
#' @export
get_sigma <- function(...) {
  UseMethod("get_sigma")
}

#' Get the covariance matrix
#' 
#' The associations in each module are taken as partial correlations, and
#' the covariance matrix is calculated from these assuming that expression
#' for gene i is the weighted average over each module using 1/sqrt(m_i) 
#' as the weight, where m_i is the number of modules containing gene i.
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @note If a matrix is provided, it is assumed to be a partial correlation matrix.
#' A warning is given in this case. To avoid the warning message, convert the
#' matrix into a network object using 'create_network_from_association_matrix()'.
#' @return A covariance matrix.
#' @export
get_sigma.default <- function(...) {
  cat("get_sigma() is defined for 'network' and 'network_module' objects.\n")
}

#' Get the covariance matrix
#' 
#' The associations in each module are taken as partial correlations, and
#' the covariance matrix is calculated from these assuming that expression
#' for gene i is the weighted average over each module using 1/sqrt(m_i) 
#' as the weight, where m_i is the number of modules containing gene i.
#' @param m A 'matrix' object.
#' @note The matrix is assumed to be a partial correlation matrix.
#' A warning is given in this case. To avoid the warning message, convert the
#' matrix into a network object using 'create_network_from_association_matrix()'.
#' @return A covariance matrix.
#' @export
get_sigma.matrix <- function(m, ...) {
  if(any(diag(m) == 0)) {
    stop("Matrix 'm' should have nonzero entries along its diagonal.")
  }
  warning("Interpreting matrix 'm' as a partial correlation matrix.")
  
  precision_matrix <- -m
  if(all(precision_matrix == 0)) {
    # If there are no connections in m, return the identify matrix.
    return(diag(1, nrow(precision_matrix)))
  }
  diag(precision_matrix) <- 1
  if(any(eigen(precision_matrix)$values < 0)) {
    warning(paste("The edge weights in the module do not correspond to a", 
                  "positive definite precision matrix."))
  }
  sigma <- solve(precision_matrix)
  
  return(sigma)
}




#' Check if an object is weighted
#' 
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return Returns FALSE if all of the connections in the 
#' object are weighted by 0s and 1s, and returns TRUE otherwise. If there
#' are no connections in the module, then this function returns TRUE.
#' @export
is_weighted <- function(...) {
  UseMethod("is_weighted")
}

#' Check if an object is weighted
#' 
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return Returns FALSE if all of the connections in the 
#' object are weighted by 0s and 1s, and returns TRUE otherwise. If there
#' are no connections in the module, then this function returns TRUE.
#' @export
is_weighted.default <- function(...) {
  cat("is_weighted() is defined for 'network' and 'network_module' objects.\n")
}

#' Check if an object is weighted
#' 
#' @param m A 'matrix' object.
#' @return Returns FALSE if all of entries in the matrix are either zero or one, 
#' and returns TRUE otherwise. However, if there are no connections, i.e. all
#' off-diagonal entries are zero, then this function returns TRUE.
#' @export
is_weighted.matrix <- function(m, ...) {
  if(!(class(m) == "matrix")) 
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a 'matrix' object."))
  if(dim(m)[1] != dim(m)[2])
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a square matrix."))
  if(max(abs(m - t(m))) > 10^-13)   
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a symmetric matrix."))
  
  m <- m[lower.tri(m)]
  
  if(any(!(m %in% c(0, 1)))) {
    return(TRUE)
  } else if(all(m == 0)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Removes the weights of all connections
#' 
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return The modified object.
#' @export
remove_weights <- function(...) {
  UseMethod("remove_weights")
}

#' Removes the weights of all connections
#' 
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return The modified object.
#' @export
remove_weights.default <- function(...) {
  cat("remove_weights() is defined for 'network', 'network_module', and 'matrix' objects.\n")
}

#' Removes the weights of all connections
#' 
#' @param m A 'matrix' object.
#' @return The modified object.
#' @export
remove_weights.matrix <- function(m, ...) {
  if(!(class(m) == "matrix")) 
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a 'matrix' object."))
  if(dim(m)[1] != dim(m)[2])
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a square matrix."))
  if(max(abs(m - t(m))) > 10^-13)   
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a symmetric matrix."))
  
  m[m != 0] <- 1
  return(m)
}



#' Get node names
#' 
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return A vector containing the node names or node indicies.
#' @export
get_node_names <- function(...) {
  UseMethod("get_node_names")
}

#' Get node names
#' 
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return A vector containing the node names or node indicies.
#' @export
get_node_names.default <- function(...) {
  cat("get_node_names() is defined for 'network' and 'network_module' objects.\n")
}

#' Get node names
#' 
#' @param m A 'matrix' object.
#' @return A vector containing the node names or node indicies.
#' @export
get_node_names.matrix <- function(m, ...) {
  node_names <- colnames(m)
  if(is.null(node_names)) {
    node_names <- 1:ncol(m)
  }
  return(node_names)
}




#' Rewire connections to a node.
#' 
#' @param ... A 'network', 'network_module', or 'matrix' object to modify. 
#' @return The modified object m.
#' @export
rewire_connections_to_node <- function(...) {
  UseMethod("rewire_connections_to_node")
}

#' Rewire connections to a node.
#' 
#' @param ... A 'network', 'network_module', or 'matrix' object to modify. 
#' @return The modified object m.
#' @export
rewire_connections_to_node.default <- function(...) {
  cat("rewire_connections_to_node() is defined for 'network', 'network_module'",
      " and 'matrix' objects.\n")
}

#' Rewire connections to a node.
#' 
#' @param m A 'matrix' object to modify. 
#' @param node The node to rewire.
#' @param rewire_prob A value between 0 and 1. Each connection to 'node' 
#' will be rewired with probability equal to 'rewire_prob'. Note, the degree of 
#' 'node' is unchanged after this operation.
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling nodes to rewire.
#' @param exponent The exponent used for weighted sampling. When exponent = 0,
#' nodes are sampled uniformly. When exponent > 0, the sampling probability
#' is based on node weights.
#' @return The modified object m.
#' @export
rewire_connections_to_node.matrix <- function(m,
                                              node,
                                              rewire_prob,
                                              weights = NULL,
                                              exponent = 0) {
  if(!(class(m) == "matrix")) 
    stop(paste0("'", deparse(substitute(m)), 
                "' is not a 'matrix' object."))
  
  nodes <- get_node_names(m)
  p <- length(nodes)
  
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'm'.
  check_adjacency_matrix(m, checklist)
  
  # Check 'node'.
  check_positive_integer(node, checklist)
  
  # Check 'rewire_prob'. (Allowed to equal 0 or 1.)
  check_in_closed_interval(rewire_prob, checklist, 0, 1)
  
  # Check 'weights'.
  if(!is.null(weights) && (length(weights) != p || !is.numeric(weights))) 
    stop(paste0("'", deparse(substitute(weights)), 
                "' is not a numeric vector of length equal to the number of",
                " nodes in m."))
  
  if(length(exponent) != 1 || !is.numeric(exponent)) 
    stop(paste0("'", deparse(substitute(exponent)), 
                "' is not positive numeric value."))
  
  report_checks(checklist)
  ##################################
  
  if(is.null(weights)) {
    weights <- rep(0, p)
  }
  
  adj_matrix <- m
  
  if(node %in% nodes) {
    node_index <- which(node == nodes)
    degree <- apply(adj_matrix, 2, sum)
    weights <- weights + degree
    
    # Rewire if the node has neighbors and is not connected to all other nodes.
    if(degree[node_index] > 0 && degree[node_index] != (p - 1)) {
      neighbors <- which(adj_matrix[, node_index] != 0)
      available <- setdiff((1:p)[-node_index], neighbors)
      
      # Determine how many connections to rewire
      n_rewire <- sum(runif(length(neighbors)) < rewire_prob)
      if(n_rewire > 0) {
        
        # Sample the connections to rewire.
        if(length(neighbors) == 1) {
          # Do not use sample() if only one neighbor.
          neighbors <- neighbors
        } else {
          neighbors <- sample(neighbors, n_rewire,
                              prob = 1 - ecdf(weights[neighbors])(weights[neighbors])^exponent + 0.001)
        }
        
        # Rewire each connection.
        for(old_neighbor in neighbors) {
          # Sample a new neighbor to wire to from the available nodes.
          if(length(available) == 1) {
            # Do not use sample() if only one available.
            new_neighbor <- available
          } else {
            new_neighbor <- sample(available, 1,
                                   prob = ecdf(weights[available])(weights[available])^exponent + 0.001)
          }
          # Rewire connection to new node.
          adj_matrix[old_neighbor, node_index] <- 0
          adj_matrix[node_index, old_neighbor] <- 0
          adj_matrix[new_neighbor, node_index] <- 1
          adj_matrix[node_index, new_neighbor] <- 1
          # The new neighbor is no longer available. Replace it with the old 
          # neighbor, which can now be selected.
          available[which(available == new_neighbor)] <- old_neighbor
          # Update weights.
          weights[new_neighbor] <- weights[new_neighbor] + 1
          weights[old_neighbor] <- weights[old_neighbor] - 1
        }
      }
      
      # Update the module structure.
      m <- adj_matrix
    }
  } else {
    warning("'", deparse(substitute(node)), "' is not a node in 'm'.",
            "Returning 'm' unmodified.")
  }
  
  return(m)
}








#' Remove connections to a node.
#' 
#' @param ... A 'network', 'network_module', or 'matrix' object to modify. 
#' @return The modified object m.
#' @export
remove_connections_to_node <- function(...) {
  UseMethod("remove_connections_to_node")
}

#' Remove connections to a node.
#' 
#' @param ... A 'network', 'network_module', or 'matrix' object to modify. 
#' @return The modified object.
#' @export
remove_connections_to_node.default <- function(...) {
  cat("remove_connections_to_node() is defined for 'network', 'network_module'",
      " and 'matrix' objects.\n")
}

#' Remove connections to a node.
#' 
#' @param adj_matrix An adjacency matrix to modify.
#' @param node The node to unwire.
#' @param remove_prob A value between 0 and 1. Each connection to 'node_index' 
#' will be removed with probability equal to 'remove_prob'.
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling neighbors to unwire from.
#' @param exponent The exponent used for weighted sampling. When exponent = 0,
#' neighboring nodes are sampled uniformly. When exponent > 0, the sampling 
#' probability is based on node weights.
#' @return The modified adjacency matrix.
#' @export
remove_connections_to_node.matrix <- function(adj_matrix,
                                              node,
                                              remove_prob,
                                              weights = NULL,
                                              exponent = 0) {
  if(!(class(adj_matrix) == "matrix")) 
    stop(paste0("'", deparse(substitute(adj_matrix)), 
                "' is not a 'matrix' object."))
  
  nodes <- get_node_names(adj_matrix)
  p <- length(nodes)
  
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'adj_matrix'.
  check_adjacency_matrix(adj_matrix, checklist)
  
  # Check 'node'.
  check_positive_integer(node, checklist)
  
  # Check 'remove_prob'. (Allowed to equal 0 or 1.)
  check_in_closed_interval(remove_prob, checklist, 0, 1)
  
  # Check 'weights'.
  if(!is.null(weights) && (length(weights) != p || !is.numeric(weights))) 
    stop(paste0("'", deparse(substitute(weights)), 
                "' is not a numeric vector of length equal to the number of",
                " nodes in 'adj_matrix'"))
  
  if(length(exponent) != 1 || !is.numeric(exponent)) 
    stop(paste0("'", deparse(substitute(exponent)), 
                "' is not positive numeric value."))
  
  report_checks(checklist)
  ##################################
  
  if(is.null(weights)) {
    weights <- rep(0, p)
  }
  
  if(node %in% nodes) {
    node_index <- which(node == nodes)
    degree <- apply(adj_matrix, 2, sum)
    weights <- weights + degree
    
    # Rewire if the node has neighbors.
    if(degree[node_index] > 0 && degree[node_index] != (p - 1)) {
      neighbors <- which(adj_matrix[, node_index] != 0)
      n_remove <- sum(runif(length(neighbors)) < remove_prob)
      
      if(n_remove > 0) {
        
        # Sample the connections to rewire.
        if(length(neighbors) == 1) {
          # Do not use sample() if only one neighbor.
          neighbors <- neighbors
        } else {
          neighbors <- sample(neighbors, n_remove,
                              prob = 1 - ecdf(weights[neighbors])(weights[neighbors])^exponent + 0.001)
        }
        
        # Remove connections
        adj_matrix[neighbors, node_index] <- 0
        adj_matrix[node_index, neighbors] <- 0
      }
    }
  } else {
    warning("'", deparse(substitute(node)), "' is not a node in 'adj_matrix'.",
            "Returning 'adj_matrix' unmodified.")
  }
  
  return(adj_matrix)
}