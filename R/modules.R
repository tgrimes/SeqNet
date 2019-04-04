library(methods)
setClass(Class = "network_module")


#' Create a module
#' 
#' @param nodes A numeric vector indicating which nodes in the network are
#' contained in this module.
#' @return A 'network_module' object.
#' @export
create_empty_module <- function(nodes) {
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check nodes
  check_positive_integer_vector(nodes, checklist)
  
  report_checks(checklist)
  ##################################

  n_nodes = length(nodes)
  
  module <- list(name = NULL,
                 nodes = nodes,
                 struct = rep(0, n_nodes * (n_nodes - 1) / 2))
  
  class(module) <- "network_module"
  
  return(module)
}



###########################################################################
#
# Functions to create a module
#
###########################################################################

#' Create a module from an adjacency matrix
#' 
#' @param nodes A numeric vector indicating which nodes in the network are
#' contained in this module.
#' @param module_name (optional) Character string specifying the name of the 
#' module. If NULL, the module will be unnamed.
#' @return A 'network_module' object.
#' @export
create_module_from_adjacency_matrix <- function(adjacency_matrix, 
                                                nodes = NULL,
                                                module_name = NULL) {
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'adjacency_matrix'
  check_adjacency_matrix(adjacency_matrix, checklist)
  
  # Check 'nodes', if provided.
  if(!is.null(nodes)) {
    check_positive_integer_vector(nodes, checklist)
    if(length(nodes) != ncol(adjacency_matrix)) 
      ArgumentCheck::addError(
        msg = "Length of argument 'nodes' must equal the number of columns of 'adjacency_matrix'.",
        argcheck = checklist
      )
  }
  
  # Check 'module_name'
  if(!is.null(module_name) && !is.character(module_name)) 
    ArgumentCheck::addError(
      msg = "Argument 'module_name' must be a character string.",
      argcheck = checklist
    )
  
  report_checks(checklist)
  ##################################
  
  if(is.null(nodes)) {
    if(!is.null(colnames(adjacency_matrix))) {
      nodes <- as.numeric(colnames(adjacency_matrix))
    } else {
      nodes <- 1:ncol(adjacency_matrix)
    }
  }
  
  module <- create_empty_module(nodes)
  module <- set_module_name(module, module_name)
  module <- set_module_struct(module, adjacency_matrix)
  
  return(module)
}




#' Create a module from an association matrix
#' 
#' The edge weights in the module will be set to the corresponding values
#' in the association matrix. Note, only the off-diagonal elements of the 
#' association matrix are used.
#' @param association matrix
#' @param nodes A numeric vector indicating which nodes in the network are
#' contained in this module.
#' @param module_name (optional) Character string specifying the name of the 
#' module. If NULL, the module will be unnamed.
#' @return A 'network_module' object.
#' @export
create_module_from_association_matrix <- function(association_matrix, 
                                                  nodes = NULL,
                                                  module_name = NULL) {
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'association_matrix'
  check_association_matrix(association_matrix, checklist)
  
  # Check 'nodes', if provided.
  if(!is.null(nodes)) {
    check_positive_integer_vector(nodes, checklist)
    if(length(nodes) != ncol(association_matrix)) 
      ArgumentCheck::addError(
        msg = "Length of argument 'nodes' must equal the number of columns of 'association_matrix'.",
        argcheck = checklist
      )
  }
  
  # Check 'module_name'
  if(!is.null(module_name) && !is.character(module_name)) 
    ArgumentCheck::addError(
      msg = "Argument 'module_name' must be a character string.",
      argcheck = checklist
    )
  
  report_checks(checklist)
  ##################################
  
  if(is.null(nodes)) {
    if(!is.null(colnames(association_matrix))) {
      nodes <- as.numeric(colnames(association_matrix))
    } else {
      nodes <- 1:ncol(association_matrix)
    }
  }
  
  adjacency_matrix <- association_matrix
  adjacency_matrix[adjacency_matrix != 0] <- 1
  
  module <- create_module_from_adjacency_matrix(adjacency_matrix, 
                                                nodes = nodes, 
                                                module_name = module_name)
  module <- set_module_weights(module, association_matrix)

  return(module)
}

#' Create a random module
#' 
#' @param nodes A numeric vector indicating which nodes in the network 
#' are contained in this module. 
#' @param module_name (optional) Character string specifying the name of the 
#' module. If NULL, the module will be unnamed.
#' @param add_weights If TRUE, weights will be generated for the module 
#' connections.
#' @param ... Additional arguments passed to 'update_module_with_random_struct()' 
#' and 'update_module_with_random_weights()'.
#' @details See ?igraph::watts.strogatz.game for details on generating
#' small-world networks.
#' @return A 'network_module' object.
#' @export
random_module <- function(nodes, 
                          module_name = NULL,
                          add_weights = FALSE,
                          ...) {
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'nodes', if provided.
  check_positive_integer_vector(nodes, checklist)
  
  report_checks(checklist)
  ##################################
  
  module <- create_empty_module(nodes)
  module <- set_module_name(module, module_name)
  module <- update_module_with_random_struct(module, ...)
  if(add_weights) {
    module <- update_module_with_random_weights(module, ...)
  }
  
  return(module)
}





###########################################################################
#
# Setter functions for network_module.
#
###########################################################################

#' Set the connections for a module
#' 
#' @param module The 'network_module' object to modify.
#' @param struct A vector or matrix of weights for each connetions. If a vector,
#' its length must equal #nodes(#nodes-1)/2, where p is the number of nodes in the module. 
#' If a matrix, it should be square with the number of columns equal to p; only
#' the upper-triangle of the matrix will be used.
#' @return The modified 'network_module' object.
#' @export
set_module_struct <- function(module, struct) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  
  p <- length(module$nodes)
  n_possible_connections <- p * (p - 1) / 2
  if(is.vector(struct)) {
    if(length(struct) != n_possible_connections) {
      stop(paste0("Length of argument 'struct' is ", length(struct), 
                  ", but must be #nodes(#nodes - 1)/2 = ", 
                  n_possible_connections, "."))
    }
    
    # Update struct.
    module$struct <- struct
  } else if(is.matrix(struct)) {
    if(ncol(struct) != nrow(struct)) {
      stop("Argument 'struct' is not a square matrix.")
    }
    if(ncol(struct) != p) {
      stop("Number of columns for argument 'struct' does not equal number of nodes in module.")
    }
    
    struct <- struct[lower.tri(struct)]
    # Update struct.
    module$struct <- struct
  } else {
    stop("Argument 'struct' must be a vector or matrix.")
  }
  
  return(module)
}

#' Set the connection weights for a module
#' 
#' @param module The 'network_module' object to modify.
#' @param weights A vector or matrix of weights for each connetions. If a vector,
#' its length must equal the number of connections in the module. If a matrix,
#' it should be square with the number of columns equal to the number of nodes 
#' in the module; only the entries in the upper-triangle that correspond to 
#' connections in the module will be used.
#' @return The modified 'network_module' object.
#' @export
set_module_weights <- function(module, weights) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  
  n_edges <- sum(module$struct)
  index_edges <- which(module$struct != 0)
  
  if(is.vector(weights)) {
    if(length(weights) != n_edges) {
      stop(paste0("Length of argument 'weights' is ", length(weights), 
                  ", but there are ", n_edges, " connections in the network."))
    }
    if(any(weights == 0)) {
      warning("Some values in argument 'weights' are zero. This will remove module connections.")
    }
    
    # Update struct with weights.
    module$struct[index_edges] <- weights
  } else if(is.matrix(weights)) {
    if(ncol(weights) != nrow(weights)) {
      stop("Argument 'weights' is not a square matrix.")
    }
    if(ncol(weights) != length(module$nodes)) {
      stop("Number of columns for argument 'weights' does not equal number of nodes in module.")
    }
    
    weights <- weights[lower.tri(weights)]
    # Update struct with weights. 
    # Only use entries in 'weights' that correspond to edges.
    module$struct[index_edges] <- weights[index_edges]
  } else {
    stop("Argument 'weights' must be a vector or matrix.")
  }
  
  return(module)
}

#' Removes the connection weights for a module
#' 
#' @param module The 'network_module' object to modify.
#' @return The modified 'network_module' object.
#' @export
remove_weights.network_module <- function(module) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  
  if(!is_weighted(module)) {
    return(module)
  }
  module$struct[which(module$struct != 0)] <- 1
  return(module)
}

#' Set the name for a module
#' 
#' @param module The 'network_module' object to modify.
#' @param module_name A character string.
#' @return The modified 'network_module' object.
#' @export
set_module_name <- function(module, module_name) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  
  if(!is.null(module_name) &&
     (is.na(module_name) || module_name == "")) {
    # If 'module_name' is NA or otherwise empty, coerce to NULL without warning.
    # This may happen if a list of modules is partially named.
    module_name <- NULL
  }
  
  if(!is.null(module_name) && !is.character(module_name)) 
    stop("Argument 'module_name' must be a character string.")
  
  module$name <- module_name
  return(module)
}

###########################################################################
#
# Getter functions for network_module.
#
###########################################################################

#' Get adjacency matrix of a module
#' 
#' @param module The 'network_module' object to get adjacency matrix for.
#' @return An adjacency matrix with entry ij == 1 if node i and j are 
#' connected, and 0 otherwise.
#' @export
get_adjacency_matrix.network_module <- function(module, ...) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  
  n_nodes <- length(module$nodes)
  adj_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  edges <- module$struct
  # Set nonzero edges to 1, in case module is weighted.
  edges[edges != 0] <- 1
  adj_matrix[lower.tri(adj_matrix)] <- edges
  adj_matrix <- adj_matrix + t(adj_matrix)
  
  colnames(adj_matrix) <- module$nodes
  return(adj_matrix)
}

#' Get association matrix of a module
#' 
#' @param module A 'network_module' object; can be either weighted or unweighted.
#' @return An association matrix with entry ij != 0 if node i and j are 
#' connected, and 0 otherwise. If the module is unweighted, then the matrix
#' will be an adjacency matrix containing ones and zeroes.
#' @export
get_association_matrix.network_module <- function(module, ...) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  
  n_nodes <- length(module$nodes)
  assoc_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  edges <- module$struct
  assoc_matrix[lower.tri(assoc_matrix)] <- edges
  assoc_matrix <- assoc_matrix + t(assoc_matrix)
  diag(assoc_matrix) <- 1
  
  # Set values near zero to exactly zero.
  assoc_matrix[abs(assoc_matrix) <= 10^-13] <- 0
  
  colnames(assoc_matrix) <- module$nodes
  return(assoc_matrix)
}

#' Get the covariance matrix of a module
#' 
#' The associations in the module are taken as partial correlations; a 
#' precision matrix is obtained using the negative partial correlations, and
#' the covariance matrix is computed by inverting the precision matrix.
#' @param module A weighted 'network_module' object.
#' @return A covariance matrix for the module
#' @export
get_sigma.network_module <- function(module, ...) {
  precision_matrix <- -get_association_matrix(module)
  if(all(precision_matrix == 0)) {
    return(diag(1, nrow(precision_matrix)))
  }
  diag(precision_matrix) <- -diag(precision_matrix)
  sigma <- solve(precision_matrix)
  return(sigma)
}


#' Get edge weights.
#' 
#' @param module The 'network_module' object to get adjacency matrix for.
#' @return A vector containing the nonzero entries of module$struct.
#' @export
get_edge_weights_from_module <- function(module) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  
  weights <- module$struct
  weights <- weights[weights != 0]
  return(weights)
}

#' Get node names of a module
#' 
#' Modules do not retain the names of each node, so the node indicies returned
#' used instead. If the network node names are known, then the vector returned
#' from this function can be used to index those node names.
#' @param module The 'network_module' object to get node names from.
#' @return A vector containing the node indicies in the module
#' @export
get_node_names.network_module <- function(module, ...) {
  return(module$nodes)
}

###########################################################################
#
# Helper functions for network_module.
#
###########################################################################

#' Generate small-world network structure for module
#' 
#' The small-world network is generated using the Watts-Strogatz method.
#' See ?igraph::watts.strogatz.game for details.
#' @param module The network_module object to modify.
#' @param lattice_neig A positive integer used for generating small-world 
#' networks for modules. If the number of nodes, p, in the module is 
#' p >= (1 + 2 * 'lattice_neig'), then
#' the module will contain (p * 'lattice_neig') edges. This can be used to 
#' the desired level of sparsity. For sparsity in [0, 1], with one indicating no 
#' edges (completely sparse), setting 
#' 'lattice_neig' = floor((1 - sparsity) * (p - 1) / 2) will constrain
#' the sparsity to above the desired level. Note, if sparsity is nearly 1, this
#' value for 'lattice_neig' may be zero, in which case p will need to be 
#' increased or sparsity will need to be lowered.
#' @param rewire_prob Used for generating small-world networks for modules; 
#' specifies the rewiring probability of each initial edges. At 0 the graph is
#' a lattice ring, and at 1 it is a random graph. 
#' @param ... Additional parameters are ignored.
#' @return An updated 'network_module' object.
#' @export
update_module_with_random_struct <- function(module, 
                                             lattice_neig = 1, 
                                             rewire_prob = 0.8, ...) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'module', if provided.
  n_nodes <- length(module$nodes)
  if(n_nodes < 2) 
    ArgumentCheck::addError(
      msg = "Argument 'module' contains 1 or fewer nodes. Cannot create struct.",
      argcheck = checklist
    )
  
  # Check 'lattice_neig'.
  if(lattice_neig <= 0) 
    ArgumentCheck::addError(
      msg = "Argument 'lattice_neig' must be greater than zero.",
      argcheck = checklist
    )
  if(lattice_neig >= n_nodes) 
    ArgumentCheck::addError(
      msg = paste0("Argument 'lattice_neig' must be less than the number",
                   " of nodes in the module (= ", n_nodes, ")."),
      argcheck = checklist
    )
  
  # Check 'reqire_prob'.
  check_in_closed_interval(rewire_prob, checklist, 0, 1)
  
  report_checks(checklist)
  ##################################
  
  p <- length(module$nodes)
  adjacency_matrix <- as.matrix(
    igraph::get.adjacency(
      igraph::watts.strogatz.game(dim = 1, 
                                  size = p, 
                                  nei = lattice_neig, 
                                  p = rewire_prob)))
  
  module <- set_module_struct(module, adjacency_matrix)
  
  return(module)
}


#' Generate small-world network structure for module
#' 
#' The small-world network is generated using the Watts-Strogatz method.
#' See ?igraph::watts.strogatz.game for details.
#' @param module The network_module object to modify.
#' @param rdist A distribution function that generates random numbers. The first
#' argument should specify the number of weights to generate. By default, 
#' weights are generated uniformly from the set (-1, -0.5)U(0.5, 1).
#' @param ... Additional parameters are ignored.
#' @return An updated 'network_module' object.
#' @export
update_module_with_random_weights <- function(module, 
                                              rdist = function(n) {
                                                runif(n, 0.5, 1) * (-1)^rbinom(n, 1, 0.5)
                                              },
                                              ...) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  
  n_weights <- sum(module$struct != 0)
  if(n_weights == 0) {
    warning("Argument 'module' contains no connections. Returning module unmodified.")
  } else {
    weights <- rdist(n_weights)
    module <- set_module_weights(module, weights)
  }

  return(module)
}



#' Rewire connections to a node from a module
#' 
#' @param node_index The index of the node in the module to remove.
#' @param module The 'network_module' object to modify.
#' @param rewire_prob The node at 'node_index' has probaility equal to 
#' 'rewire_prob' to be connected with any other node in the module.
#' @return The modified 'network_module' object.
#' @export
rewire_connections_to_node_in_module <- function(node,
                                                 module,
                                                 rewire_prob = 0.2) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'node'.
  check_positive_integer(node, checklist)
  
  # Check 'rewire_prob'. (Allowed to equal 0 or 1.)
  check_in_closed_interval(rewire_prob, checklist, 0, 1)
  
  report_checks(checklist)
  ##################################
  
  if(node %in% module$nodes) {
    p <- length(module$nodes)
    node_index <- which(node == module$nodes)
    adj_matrix <- get_adjacency_matrix(module)
    # Generate connections to other nodes with probability = 'rewire_prob'.
    new_connections <- rbinom(p, 1, rewire_prob)
    # No connection to self.
    new_connections[node_index] <- 0
    # Update adjacency matrix.
    adj_matrix[node_index, ] <- new_connections
    adj_matrix[, node_index] <- new_connections
    # Update module with new connections.
    module <- set_module_struct(module, adj_matrix)
  } else {
    warning("'", deparse(substitute(node)), "' is not a node in 'module'.",
            "Returning module unmodified.")
  }
  
  return(module)
}

#' Subset module onto largest component
#' 
#' Any nodes that are disconnected from the largest component in the module
#' are removed.
#' @param module The 'network_module' object to modify.
#' @return The modified 'network_module' object.
#' @export
remove_small_components_from_module <- function(module) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  
  adj_matrix <- get_adjacency_matrix(module)
  g <- igraph::graph_from_adjacency_matrix(adj_matrix,
                                           mode = "undirected")
  comp <- igraph::components(g)
  if(comp$no > 1) {
    # If there are more than one component in the graph, subset onto the largest
    # one. Use association matrix to preserve weights, if they are present.
    largest_component <- which(comp$csize == max(comp$csize))[1]
    keep_nodes <- which(comp$membership == largest_component)
    assoc_matrix <- get_association_matrix(module)
    assoc_matrix <- assoc_matrix[keep_nodes, keep_nodes]
    module <- create_module_from_association_matrix(assoc_matrix)
  } 
  
  return(module)
}

###########################################################################
#
# Utility functions for network_module.
#
###########################################################################

#' Print function for 'network_module' object.
#' 
#' @param module A 'network_module' object.
#' @param ... Additional arguments are ignored.
#' @return Prints a summary of the module.
#' @export
print.network_module <- function(module, ...) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  module_message <- 
    paste0("Module: ", 
           ifelse(is.null(module$name), 
                  "-unnamed-", 
                  module$name),
           "\n")
  
  n_nodes <- length(module$nodes)
  node_message <- 
    paste0("Contains nodes:",
           ifelse(n_nodes > 50, 
                  paste(paste(module$nodes[1:50], collapse = ", "), "..."),
                  paste(module$nodes, collapse = ", ")),
           "\n",
           "# nodes: ", n_nodes, "\n")
  
  n_edges <- sum(module$struct != 0)
  edge_message <- 
    paste0("# edges: ", n_edges, "\n")
  
  weighted <- !all(module$struct %in% c(0, 1))
  weight_message <- 
    ifelse(weighted,
           paste("Connection are weighted"),
           paste("Connections are unweighted"))
  
  message <- paste0(module_message, node_message, edge_message, weight_message)
  cat(message)
}

#' Check if a module is weighted
#' 
#' @param module The 'network_module' object to check.
#' @return A boolean value that is FALSE if all of the connections in the 
#' module are weighted by 0s and 1s, and it returns TRUE otherwise. If there
#' are no connections in the module, then this function returns TRUE.
#' @export
is_weighted.network_module <- function(module, ...) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  if(all(module$struct == 0)) {
    return(TRUE)
  }
  if(any(!(module$struct %in% c(0, 1)))) {
    return(TRUE)
  } 
  return(FALSE)
}