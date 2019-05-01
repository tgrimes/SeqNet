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
  
  nodes <- sort(nodes)
  
  module <- list(name = NULL,
                 nodes = nodes,
                 edges = NULL)
  
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
#' The edges in the module will be set to the edges in the adjacency matrix. 
#' The edges are undirected, and only the lower triangle of the
#' matrix is considered. See ?set_module_edges for more details. 
#' @param adjacency_matrix The adjacency matrix used to create the module.
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
  module <- set_module_edges(module, adjacency_matrix)
  
  return(module)
}




#' Create a module from an association matrix
#' 
#' The edge weights in the module will be set to the corresponding values
#' in the association matrix. The edges are undirected, and only the lower 
#' triangle of the matrix is considered. See ?set_module_weights for more details. 
#' @param association_matrix The association matrix used to create the module.
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
#' @param ... Additional arguments passed to 'update_module_with_random_edges()' 
#' and 'update_module_with_random_weights()'.
#' @details See ?igraph::watts.strogatz.game for details on generating
#' small-world networks.
#' @return A 'network_module' object.
#' @export
random_module <- function(nodes, 
                          module_name = NULL,
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
  module <- update_module_with_random_edges(module, ...)
  
  return(module)
}





###########################################################################
#
# Setter functions for network_module.
#
###########################################################################

#' Set the edges in a module
#' 
#' @param module The 'network_module' object to modify.
#' @param edges A matrix used to indicate the edges in the module. If the matrix
#' is square and contains the same number of rows and columns as nodes in 
#' the module, then it is assumed to be an adjacency matrix and the nonzero 
#' lower-triangle values of the matrix are used to indicate edges in the module.
#' If the matrix is not square, the first two columns are assumed to be an
#' edge list.
#' @return The modified 'network_module' object.
#' @export
set_module_edges <- function(module, edges) {
  if(!(class(module) == "network_module")) 
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  # If there is only one node in the module, return it unmodified.
  if(length(module$nodes) == 1) {
    return(module)
  }
  
  if(is.null(edges)) {
    module$edges <- NULL
  } else if(is.matrix(edges)) {
    # If a square matrix is provided, 
    if(nrow(edges) == ncol(edges)) {
      if(nrow(edges) != length(module$nodes)) {
        stop(paste0("Argument 'edges' is a square matrix, but the number of", 
                    "columns does not match number of nodes in the module."))
      }
      # Set any nonzero values to 1. Use only the lower-triangle entries.
      edges[edges != 0] <- 1
      edges[upper.tri(edges)] <- 0
      diag(edges) <- 0
      edges <- edges + t(edges)
      colnames(edges) <- 1:ncol(edges) # Use indicies in edge list (edges).
      edges <- igraph::graph.adjacency(edges, mode = "undirected")
      edges <- apply(igraph::get.edgelist(edges), 2, as.numeric)
    }
    if(nrow(edges) < 2) {
      stop("Argument 'edges' must be a matrix with at least two columns.")
    }
    # Update edges.
    module$edges <- edges
  } else {
    stop("Argument 'edges' must be a matrix.")
  }
  
  return(module)
}

#' Set the connection weights for a module
#' 
#' @param module The 'network_module' object to modify.
#' @param weights A vector or matrix of weights for each connetions. If a vector,
#' its length must equal the number of connections in the module. If a matrix,
#' it should be square with the number of columns equal to the number of nodes 
#' in the module; only the entries in the lower triangle that correspond to 
#' connections in the module will be used.
#' @return The modified 'network_module' object.
#' @export
set_module_weights <- function(module, weights) {
  if(!(class(module) == "network_module")) 
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  if(is.null(module$edges) && !is.null(weights)) {
    warning("Argument 'weights' is not NULL, but the module edges is.")
    # Nothing to do, so return module unchanged.
    return(module)
  }
  
  n_edges <- nrow(module$edges)
  
  if(is.vector(weights)) {
    if(length(weights) != n_edges) {
      stop(paste0("Length of argument 'weights' is ", length(weights), 
                  ", but there are ", n_edges, " connections in the network."))
    }
  } else if(is.matrix(weights)) {
    if(ncol(weights) != nrow(weights)) {
      stop("Argument 'weights' is not a square matrix.")
    }
    if(ncol(weights) != length(module$nodes)) {
      stop(paste0("Argument 'weights' is a square matrix, but the number of", 
                  "columns does not match number of nodes in the module."))
    }
    # Get values in lower-triangle that correspond to connections in the module.
    weights <- weights[module$edges[, 2:1]]
  } else {
    stop("Argument 'weights' must be a vector or matrix.")
  }
  
  # Update edges with weights.
  module$edges <- cbind(module$edges, weights)
  
  return(module)
}

#' Removes the connection weights for a module
#' 
#' @param module The 'network_module' object to modify.
#' @return The modified 'network_module' object.
#' @export
remove_weights.network_module <- function(module) {
  if(!(class(module) == "network_module")) 
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  if(is.null(module$edges)) {
    return(module)
  }
  if(!is_weighted(module)) {
    return(module)
  }
  
  # Reset edges with only first two columns.
  module$edges <- module$edges[, 1:2]
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
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
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
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  n_nodes <- length(module$nodes)
  adj_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  if(!is.null(module$edges)) {
    adj_matrix[module$edges[, 1:2]] <- 1
    adj_matrix[module$edges[, 2:1]] <- 1
  }
  
  colnames(adj_matrix) <- module$nodes
  return(adj_matrix)
}

#' Get association matrix of a module
#' 
#' @param module A weighted 'network_module' object.
#' @return An association matrix with entry ij != 0 if node i and j are 
#' connected, and 0 otherwise. The diagonal of the association matrix is set to 0.
#' @export
get_association_matrix.network_module <- function(module, ...) {
  if(!(class(module) == "network_module")) 
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  if(!is_weighted(module)) 
    stop("'", deparse(substitute(module)), "' is not weighted.")
  
  n_nodes <- length(module$nodes)
  assoc_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  if(!is.null(module$edges)) {
    # Use weights if available, otherwise set associations to 1.
    if(ncol(module$edges) >= 3) {
      weights <- module$edges[, 3]
    } else {
      weights <- rep(1, nrow(module$edges))
    }
    assoc_matrix[module$edges[, 1:2]] <- weights
    assoc_matrix[module$edges[, 2:1]] <- weights
  }
  
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
    # If there are no connections in the module, return the identify matrix.
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


#' Get edge weights.
#' 
#' @param module The 'network_module' object to get edge weights for.
#' @return A vector containing the weights of each edge. If the edges are 
#' unweighted, then a vector of 1's is returned. If there are no edges, in the
#' module, then NULL is returned.
#' @export
get_edge_weights_from_module <- function(module) {
  if(!(class(module) == "network_module")) 
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  if(is.null(module$edges)) {
    return(NULL)
  }
  if(ncol(module$edges) >= 3) {
    weights <- module$edges[, 3]
  } else {
    weights <- rep(1, nrow(module$edges))
  }
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
update_module_with_random_edges <- function(module, 
                                            ...) {
  if(!(class(module) == "network_module")) 
    stop("'", deparse(substitute(module)), 
         "' is not a 'network_module' object.")
  
  adjacency_matrix <- random_module_structure(length(module$nodes), ...)
  adjacency_matrix <- connect_module_structure(adjacency_matrix, ...)
  module <- set_module_edges(module, adjacency_matrix)
  return(module)
}

random_module_structure <- function(size, 
                                    p_rewire = 0.5,
                                    p_remove = 0.5,
                                    deg_ex = NULL,
                                    neig_size = 3,
                                    neig_size_fn = NULL,
                                    exponent = 100,
                                    ...) {
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  if(exponent < 0) 
    ArgumentCheck::addError(
      msg = "Argument 'exponent' must be positive.",
      argcheck = checklist
    )
  
  if(!is.null(neig_size_fn)) {
    neig_size <- neig_size_fn(size)
  }
  
  if(neig_size < 0) 
    ArgumentCheck::addError(
      msg = "'neig_size' must be positive..",
      argcheck = checklist
    )
  
  if(p_rewire < 0 || p_rewire > 1) 
    ArgumentCheck::addError(
      msg = "Argument 'p_rewire' must be between 0 and 1.",
      argcheck = checklist
    )
  if(p_remove < 0 || p_remove > 1) 
    ArgumentCheck::addError(
      msg = "Argument 'p_remove' must be between 0 and 1.",
      argcheck = checklist
    )
  
  report_checks(checklist)
  ##################################
  
  nodes <- 1:size
  
  adj <- as.matrix(igraph::as_adjacency_matrix(
    igraph::make_lattice(length = size, dim = 1, nei = neig_size,
                         circular = TRUE)),
    size, size)
  
  deg <- apply(adj, 2, sum)
  for(i in sample(nodes)) {
    adj <- rewire_connections_to_node(adj, i, p_rewire, deg_ex, exponent)
    adj <- remove_connections_to_node(adj, i, p_remove, deg_ex, exponent)
  }
  
  return(adj)
}

connect_module_structure <- function(adj,
                                     deg_ex = NULL,
                                     ...) {
  deg <- apply(adj, 2, sum)
  nodes <- 1:length(deg)
  if(is.null(deg_ex)) {
    deg_ex <- rep(0, ncol(adj))
  }
  
  # If there multiple components, combine into one.
  components <- igraph::components(
    igraph::graph_from_adjacency_matrix(adj, mode = "undirected"))
  if(components$no > 1) {
    main_component <- which(components$membership == 1)
    for(i in 2:components$no) {
      sub_component <- which(components$membership == i)
      # Choose a representative for the subcomponent.
      if(length(sub_component) == 1) {
        index_rep <- sub_component
      } else {
        index_rep <- sample(sub_component, 1) 
      }
      # Choose a representative for the main component.
      if(length(main_component) == 1) {
        index_main <- main_component
      } else {
        index_main <- sample(main_component, 1, 
                             prob = ecdf(deg[main_component] + 
                                           deg_ex[main_component])(deg[main_component] + 
                                                                     deg_ex[main_component])^10 + 0.001)
      }
      
      adj[index_main, index_rep] <- 1
      adj[index_rep, index_main] <- 1
      deg[c(index_rep, index_main)] <- deg[c(index_rep, index_main)] + 1
      
      main_component <- union(main_component, sub_component)
    }
  }
  
  return(adj)
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
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  if(is.null(module$edges)) {
    warning("Argument 'module' contains no connections. Returning module unmodified.")
  }else {
    weights <- rdist(nrow(module$edges))
    module <- set_module_weights(module, weights)
  }

  return(module)
}



#' Rewire connections to a node.
#' 
#' @param module A 'network_module' object to modify. 
#' @param node The node to rewire.
#' @param rewire_prob A value between 0 and 1. Each connection to 'node' 
#' will be rewired with probability equal to 'rewire_prob'. Note, the degree of 
#' 'node' is unchanged after this operation.
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling nodes to rewire.
#' @param exponent The exponent used for weighted sampling. When exponent = 0,
#' nodes are sampled uniformly. When exponent > 0, the sampling probability
#' is based on node weights.
#' @return The modified module.
#' @export
rewire_connections_to_node.network_module <- function(module,
                                                      node,
                                                      rewire_prob,
                                                      weights = NULL,
                                                      exponent = 0) {
  if(!(class(module) == "network_module")) 
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  adj_matrix <- get_adjacency_matrix(module)
  adj_matrix <- rewire_connections_to_node(adj_matrix, node, 
                                           rewire_prob, weights, exponent)
  module <- set_module_edges(module, adj_matrix)
  
  return(module)
}


#' Remove connections to a node.
#' 
#' @param module A 'network_module' object to modify. 
#' @param node The node to unwire.
#' @param remove_prob A value between 0 and 1. Each connection to 'node_index' 
#' will be removed with probability equal to 'remove_prob'.
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling neighbors to unwire from.
#' @param exponent The exponent used for weighted sampling. When exponent = 0,
#' neighboring nodes are sampled uniformly. When exponent > 0, the sampling 
#' probability is based on node weights.
#' @return The modified module.
#' @export
remove_connections_to_node.network_module <- function(module,
                                                      node,
                                                      remove_prob,
                                                      weights = NULL,
                                                      exponent = 0) {
  if(!(class(module) == "network_module")) 
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  adj_matrix <- get_adjacency_matrix(module)
  adj_matrix <- remove_connections_to_node(adj_matrix, node, 
                                           rewire_prob, weights, exponent)
  module <- set_module_edges(module, adj_matrix)
  
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
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  if(is.null(module$edges)) {
    return(module)
  }
  
  g <- igraph::graph_from_edgelist(module$edges, directed = FALSE)
  comp <- igraph::components(g)
  if(comp$no > 1) {
    # If there are more than one component in the graph, subset onto the largest
    # one. Use association matrix to preserve weights, if they are present.
    largest_component <- which(comp$csize == max(comp$csize))[1]
    keep_nodes <- which(comp$membership == largest_component)
    adj_matrix <- get_adjacency_matrix(module)
    adj_matrix <- adj_matrix[keep_nodes, keep_nodes]
    module <- create_module_from_adjacency_matrix(adj_matrix)
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
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  vals <- get_network_characteristics(module)
  message <- paste0(ifelse(is_weighted(module), "A weighted", "An unweighted"), 
                    " module containing ", vals$p, " nodes and ", 
                    vals$n_edges, " edges.\n",
                    "Contains nodes: ",
                    ifelse(vals$p > 50, 
                           paste(paste(module$nodes[1:50], collapse = ", "), "..."),
                           paste(module$nodes, collapse = ", ")), "\n")
  
  cat(message)
  print(round(unlist(vals[-c(1, 2, 6)]), 3))
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
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  
  if(is.null(module$edges)) {
    return(TRUE)
  }
  if(ncol(module$edges) >= 3) {
    return(TRUE)
  } 
  return(FALSE)
}