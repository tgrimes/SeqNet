#' Create a module
#' 
#' @param nodes A numeric vector indicating which nodes in the network are
#' contained in this module.
#' @return A 'network_module' object.
#' @export
#' @examples 
#' module <- create_empty_module(1:10)
#' plot(module) # A module with no edges.
create_empty_module <- function(nodes) {
  if(any(nodes <= 0)) 
    stop("Argument 'nodes' must contain positive integers.")
  
  if(length(nodes) != length(unique(nodes))) {
    stop("Duplicate nodes are present in the module.")
  }
  
  if(is.unsorted(nodes)) {
    nodes <- sort(nodes)
  }
  
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
#' @param run_checks If TRUE, then the 'adjacency_matrix' argument is checked.
#' @return A 'network_module' object.
#' @export
#' @examples
#' nw <- random_network(10)
#' nw <- gen_partial_correlations(nw)
#' adj_mat <- get_adjacency_matrix(nw)
#' create_module_from_adjacency_matrix(adj_mat)
create_module_from_adjacency_matrix <- function(adjacency_matrix, 
                                                nodes = NULL,
                                                module_name = NULL,
                                                run_checks = TRUE) {
  if(run_checks && !is_adjacency_cpp(adjacency_matrix)) {
    if(!all(diag(adjacency_matrix) == 0))
      stop("Argument 'adjacency_matrix' has nonzero diagonal values.")
    if(!all(abs(adjacency_matrix - t(adjacency_matrix)) < 10^-13))
      stop("Argument 'adjacency matrix' is not symmetric.")
    if(!all(adjacency_matrix %in% c(0, 1)))
      stop("Argument 'adjacency_matrix' contains values that are not 0 or 1.")
    stop("Argument 'adjacency_matrix' is not an adjacency matrix.")
  }
    
  if(run_checks && !is.null(nodes)) {
    if(any(nodes <= 0)) {
      stop("Argument 'nodes' must contain positive integers.")
    } else if(length(nodes) != ncol(adjacency_matrix)) {
      stop( "Length of argument 'nodes' must equal the number of columns of 'adjacency_matrix'.")
    }
  }
  if(run_checks && !is.null(module_name) && !is.character(module_name))
    stop("Argument 'module_name' must be a character string.")
  
  
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
#' @examples
#' nw <- random_network(10)
#' nw <- gen_partial_correlations(nw)
#' assoc_mat <- get_association_matrix(nw)
#' create_module_from_association_matrix(assoc_mat)
create_module_from_association_matrix <- function(association_matrix, 
                                                  nodes = NULL,
                                                  module_name = NULL) {
  if(!is.null(nodes)) {
    if(any(nodes <= 0)) {
      stop("Argument 'nodes' must contain positive integers.")
    } else if(length(nodes) != ncol(association_matrix)) {
      stop( "Length of argument 'nodes' must equal the number of columns of 'association_matrix'.")
    }
  }
  if(!is.null(module_name) && !is.character(module_name))
    stop("Argument 'module_name' must be a character string.")
  
  
  if(is.null(nodes)) {
    if(!is.null(colnames(association_matrix))) {
      nodes <- as.numeric(colnames(association_matrix))
    } else {
      nodes <- 1:ncol(association_matrix)
    }
  }
  
  # Adjust association matrix to guaruntee the precision matrix will be pos. def.
  if(!is_symmetric_cpp(association_matrix)) {
    warning(paste("The association weights is not symmetric. Using lower triangle."))
    association_matrix[upper.tri(association_matrix)] <- 0
    association_matrix <- (association_matrix + t(association_matrix))
    # Diagonal is doubled, but will be set to 1 later regardless.
  }
  
  # Convert to precision matrix when checking for positive definiteness. 
  # Use negative off-diagons, and set diagon to 1.
  association_matrix <- -association_matrix
  diag(association_matrix) <- 1
  
  make_adjustments <- FALSE
  if(any(abs(association_matrix[lower.tri(association_matrix)]) >= 1)) {
    warning(paste("The association weights are not in (-1, 1). Adjusting..."))
    make_adjustments <- TRUE
  } else if(!is_PD(association_matrix)) {
    warning(paste("The association matrix does not correspond to a", 
                  "positive definite precision matrix. Adjusting..."))
    make_adjustments <- TRUE
  }
  
  if(make_adjustments) {
    diag(association_matrix) <- 0
    lambda <- as.numeric(eigen(association_matrix, symmetric = TRUE, 
                               only.values = TRUE)$values)
    adjustment <- (max(lambda) * 10^-2.5 - min(lambda))
    diag(association_matrix) <- adjustment
    association_matrix <- cov2cor(association_matrix)
  }
  
  # Convert back to partial correlations. Set diagonal to zero.
  association_matrix <- -association_matrix
  diag(association_matrix) <- 0
  
  adjacency_matrix <- (association_matrix != 0) * 1
  
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
#' @param ... Additional arguments passed to 'random_module_structure()'.
#' @return A 'network_module' object.
#' @export
#' @examples 
#' module <- random_module(1:10)
random_module <- function(nodes, 
                          module_name = NULL,
                          ...) {
  if(any(nodes <= 0)) 
    stop("Argument 'nodes' must contain positive integers.")
  if(!is.null(module_name) && !is.character(module_name))
    stop("Argument 'module_name' must be a character string.")
  if(length(nodes) == 1) {
    warning("'", deparse(substitute(nodes)), "' has length one. ",
            "Expected a vector of node indicies. This module will ",
            "only contain one node.")
  }
  
  module <- create_empty_module(nodes)
  module <- set_module_name(module, module_name)
  
  adjacency_matrix <- random_module_structure(length(nodes), ...)
  module <- set_module_edges(module, adjacency_matrix)
  
  return(module)
}





###########################################################################
#
# Setter functions for network_module.
#
###########################################################################

#' Internal function used to set the edges in a module
#' 
#' @param module The 'network_module' object to modify.
#' @param edges A matrix used to indicate the edges in the module. If the matrix
#' is square and contains the same number of rows and columns as nodes in 
#' the module, then it is assumed to be an adjacency matrix and the nonzero 
#' lower-triangle values of the matrix are used to indicate edges in the module.
#' If the matrix is not square, the first two columns are assumed to be an
#' edge list.
#' @return The modified 'network_module' object.
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
      # Then edges is interpreted as an association matrix.
      if(nrow(edges) != length(module$nodes)) {
        stop(paste0("Argument 'edges' is a square matrix, but the number of", 
                    "columns does not match number of nodes in the module."))
      }
      edges <- edges_from_adjacency_cpp(edges)
      if(nrow(edges) == 0) {
        #If there are no edges, set edges to NULL.
        edges <- NULL
      }
    } else if(ncol(edges) < 2) {
      stop("Argument 'edges' must be a matrix with at least two columns.")
    }
    # Update edges.
    module$edges <- edges
  } else {
    stop("Argument 'edges' must be a matrix.")
  }
  
  return(module)
}

#' Internal function to set the connection weights for a module
#' 
#' @param module The 'network_module' object to modify.
#' @param weights A vector or matrix of weights for each connetions. If a vector,
#' its length must equal the number of connections in the module. If a matrix,
#' it should be square with the number of columns equal to the number of nodes 
#' in the module; only the entries in the lower triangle that correspond to 
#' connections in the module will be used.
#' @return The modified 'network_module' object.
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
    if(length(weights) != n_edges)
      stop("Length of argument 'weights' is ", length(weights), 
           ", but there are ", n_edges, " connections in the network.")
  } else if(is.matrix(weights)) {
    if(ncol(weights) != nrow(weights))
      stop("Argument 'weights' is not a square matrix.")
    if(ncol(weights) != length(module$nodes))
      stop("Argument 'weights' is a square matrix, but the number of", 
           "columns does not match number of nodes in the module.")
    # Get values in lower-triangle that correspond to connections in the module.
    if(n_edges == 1) {
      # If only 1 edge, index the weight matrix idrectly.
      weights <- weights[module$edges[1, 2], module$edges[1, 1]]
    } else {
      # Otherwise, index using the 2-column matrix of indicies.
      weights <- weights[module$edges[, 2:1]]
    }
  } else {
    stop("Argument 'weights' must be a vector or matrix.")
  }
  
  # Update edges with weights.
  if(ncol(module$edges) == 2) {
    # If module is currently unweighted, add a new column for weights.
    module$edges <- cbind(module$edges, weights = unname(weights))
  } else {
    # Otherwise, update the weights.
    module$edges[, 3] <- unname(weights)
  }
  
  return(module)
}

#' @inherit remove_weights
#' @export
remove_weights.network_module <- function(x, ...) {
  if(!(class(x) == "network_module")) 
    stop("'", deparse(substitute(x)), "' is not a 'network_module' object.")
  if(is.null(x$edges))
    return(x)
  if(!is_weighted(x))
    return(x)
  
  # Reset edges with only first two columns.
  if(nrow(x$edges) == 1) {
    x$edges <- matrix(x$edges[, 1:2], nrow = 1)
  } else {
    x$edges <- x$edges[, 1:2]
  }
  
  return(x)
}

#' Set the name for a module
#' 
#' @param module The 'network_module' object to modify.
#' @param module_name A character string.
#' @return The modified 'network_module' object.
#' @export
#' @examples 
#' nw <- random_network(10)
#' nw <- gen_partial_correlations(nw)
#' module <- nw$modules[[1]]
#' named_module <- set_module_name(module, "new name")
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

#' @inherit get_adjacency_matrix
#' @export
get_adjacency_matrix.network_module <- function(x, ...) {
  if(!(class(x) == "network_module")) 
    stop("'", deparse(substitute(x)), "' is not a 'network_module' object.")
  
  n_nodes <- length(x$nodes)
  adj_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  if(!is.null(x$edges)) {
    if(nrow(x$edges) == 1) {
      adj_matrix[x$edges[, 1], x$edges[, 2]] <- 1
      adj_matrix[x$edges[, 2], x$edges[, 1]] <- 1
    } else {
      adj_matrix[x$edges[, 1:2]] <- 1
      adj_matrix[x$edges[, 2:1]] <- 1
    }
  }
  
  colnames(adj_matrix) <- x$nodes
  return(adj_matrix)
}

#' @inherit get_association_matrix
#' @export
get_association_matrix.network_module <- function(x, tol = 10^-13, ...) {
  if(!(class(x) == "network_module")) 
    stop("'", deparse(substitute(x)), "' is not a 'network_module' object.")
  if(!is_weighted(x)) 
    stop("'", deparse(substitute(x)), "' is not weighted.")
  
  n_nodes <- length(x$nodes)
  assoc_matrix <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  if(!is.null(x$edges)) {
    # Use weights if available, otherwise set associations to 1.
    if(ncol(x$edges) >= 3) {
      weights <- x$edges[, 3]
    } else {
      weights <- rep(1, nrow(x$edges))
    }
    if(nrow(x$edges) == 1) {
      assoc_matrix[x$edges[, 1], x$edges[, 2]] <- weights
      assoc_matrix[x$edges[, 2], x$edges[, 1]] <- weights
    } else {
      assoc_matrix[x$edges[, 1:2]] <- weights
      assoc_matrix[x$edges[, 2:1]] <- weights
    }
  }
  
  assoc_matrix[abs(assoc_matrix) < tol] <- 0 # Set entries near zero to zero.
  colnames(assoc_matrix) <- x$nodes
  return(assoc_matrix)
}

#' @inherit get_sigma
#' @export
get_sigma.network_module <- function(x, ...) {
  # Set diagonal to 1 and flip sign of all off-diagonal entries.
  precision_matrix <- -get_association_matrix(x, tol = 0)
  diag(precision_matrix) <- 1
  if(!is_PD(precision_matrix)) {
    stop(paste("The edge weights in the module do not correspond to a", 
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
#' @examples
#' nw <- random_network(10)
#' nw <- gen_partial_correlations(nw)
#' module <- nw$modules[[1]]
#' get_edge_weights_from_module(module)
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

#' @inherit get_node_names
#' @export
get_node_names.network_module <- function(x, ...) {
  return(x$nodes)
}

###########################################################################
#
# Helper functions for network_module.
#
###########################################################################

#' Create a random network structure for a module
#' 
#' A single, connected graph is created. The graph is initialized as a ring 
#' lattice, and edges are randomly rewired and/or removed. The procedure
#' is similar to the Watts-Strogatz method, but the sampling of edges to 
#' modify can be based on the degree of each node.
#' @param size The number of nodes to include in the graph.
#' @param prob_rewire The probability of rewiring an edge.
#' @param prob_remove The probability of removing an edge.
#' @param weights (Optional) Weights used for sampling nodes. See 
#' ?rewire_connections_to_node and ?remove_connections_to_node for details.
#' @param neig_size The neighborhood size within which the nodes of the 
#' ring lattice are connected. The initial degree of each node is 2 * 'neig_size',
#' so long as 'size' >= (1 + 2 * 'neig_size') 
#' @param alpha A positive value used to parameterize the Beta distribution.
#' @param beta  A positive value used to parameterize the Beta distribution. 
#' @param epsilon A small constant added to the sampling probability of each node.
#' @param ... Additional arguments are ignored.
#' @return An adjacency matrix representing the network structure.
#' @export
#' @examples 
#' # Create a random module structure (an adjacency matrix) for 10 nodes.
#' adj_mat <- random_module_structure(10)
#' # A network object can be created using this structure.
#' module <- create_module_from_adjacency_matrix(adj_mat)
#' nw <- create_network_from_modules(10, module)
random_module_structure <- function(size, 
                                    prob_rewire = 1,
                                    prob_remove = 0.5,
                                    weights = NULL,
                                    neig_size = 3,
                                    alpha = 100,
                                    beta = 1,
                                    epsilon = 10^-5,
                                    ...) {
  if(size < 1)
    stop("Argument 'size' must be positive.")
  if(prob_rewire < 0 || prob_rewire > 1) 
    stop("Argument 'prob_rewire' must be between 0 and 1.")
  if(prob_remove < 0 || prob_remove > 1) 
    stop("Argument 'prob_remove' must be between 0 and 1.")
  if(neig_size < 0) 
    stop("Argument 'neig_size' must be positive.")
  if(length(size) != 1) 
    stop("Argument 'size' must be a single numeric value.")
  
  nodes <- 1:size
  
  adj <- ring_lattice_cpp(size, neig_size)
  
  # Go through each node, in random order, and rewire its edges.
  for(i in sample(nodes)) {
    adj <- rewire_connections_to_node(adj, i, prob_rewire, weights, 
                                      alpha, beta, epsilon, run_checks = FALSE, ...)
  }
  # Remove edges from the network.
  adj <- remove_connections(adj, prob_remove, run_checks = FALSE, ...)
  
  # Connect any disconnected components in the module.
  adj <- connect_module_structure(adj, weights, alpha, beta, epsilon)
  
  return(adj)
}

#' Connect disconnected components in an adjacency matrix
#' 
#' @param adj An adjacency matrix to modify.
#' @param weights (Optional) weights used for sampling nodes.
#' @param alpha A positive value used to parameterize the Beta distribution.
#' @param beta  A positive value used to parameterize the Beta distribution.
#' @param epsilon A small constant added to the sampling probability of each node.
#' @return A modified adjacency matrix
#' @note When connecting two components, a node is sampled from each with
#' probability proportional to ecdf(weights)(weights)^eta + epsilon,
#' where 'weights' are subset to only those nodes in the corresponding component.
#' When 'eta' = 0, this results in uniform sampling. When 'eta' > 0, 
#' nodes having larger 'weight' are more likely to be selected, where 'weight'
#' is equal to 'weights' + degree. (If Arugment 'weights' is NULL, then 'weight'
#' is simply the node degree).
#' @export
#' @examples
#' # This function is used in `random_module_structure()` to reconnect any
#' # disconnected components. To demonstrate, we'll create a random structure,
#' # remove connections to one of the nodes (that node will then be a disconnected
#' # component), and use `connect_module_structure()` to reconnect it back to
#' # the main component.
#' adj <- random_module_structure(10)
#' adj <- remove_connections_to_node(adj, 1, prob_remove = 1)
#' # Note that there are now two components in the network:
#' components_in_adjacency(adj) 
#' g <- plot_network(adj)
#' # After connecting, the network contains one component.
#' adj <- connect_module_structure(adj)
#' components_in_adjacency(adj) 
#' plot_network(adj, g)
connect_module_structure <- function(adj,
                                     weights = NULL,
                                     alpha = 100,
                                     beta = 1,
                                     epsilon = 10^-5) {
  nodes <- 1:ncol(adj)
  if(is.null(weights)) {
    weights <- rep(0, ncol(adj))
  }
  weights <- colSums(adj) + weights
  
  # Get connected components.
  components <- components_in_adjacency(adj) 
  csize <- unname(table(components))
  main_component <- which.max(csize)
  
  # If there multiple components, combine into one.
  if(length(csize) > 1) {
    nodes_in_main_component <- which(components == main_component)
    for(i in (1:length(csize))[-main_component]) {
      nodes_in_sub_component <- which(components == i)
      # Choose a representative for the subcomponent.
      if(length(nodes_in_sub_component) == 1) {
        index_rep <- nodes_in_sub_component
      } else {
        # index_rep <- sample(nodes_in_sub_component, 1, prob = weights[nodes_in_sub_component]^eta + epsilon)
        # index_rep <- sample(nodes_in_sub_component, 1,
        #                     prob = ecdf_cpp(weights[nodes_in_sub_component])^eta + epsilon)
        index_rep <- sample(nodes_in_sub_component, 1,
                            prob = pbeta(ecdf_cpp(weights[nodes_in_sub_component]), alpha, beta) + epsilon)
      }
      # Choose a representative for the main component.
      if(length(nodes_in_main_component) == 1) {
        index_main <- nodes_in_main_component
      } else {
        # index_main <- sample(nodes_in_main_component, 1, prob = weights[nodes_in_main_component]^eta + epsilon)
        # index_main <- sample(nodes_in_main_component, 1,
        #                      prob = ecdf_cpp(weights[nodes_in_main_component])^eta + epsilon)
        index_main <- sample(nodes_in_main_component, 1,
                             prob = pbeta(ecdf_cpp(weights[nodes_in_main_component]), alpha, beta) + epsilon)
      }
      
      adj[index_main, index_rep] <- 1
      adj[index_rep, index_main] <- 1
      weights[c(index_rep, index_main)] <- weights[c(index_rep, index_main)] + 1
      
      nodes_in_main_component <- union(nodes_in_main_component, nodes_in_sub_component)
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
#' @examples 
#' # Create a random module. 
#' module <- random_module(1:10)
#' is_weighted(module)
#' # Add a random weight to each connection.
#' module <- update_module_with_random_weights(module)
#' is_weighted(module)
update_module_with_random_weights <- function(module, 
                                              rdist = function(n) {
                                                runif(n, 0.5, 1) * (-1)^rbinom(n, 1, 0.5)
                                              },
                                              ...) {
  if(!(class(module) == "network_module")) 
    stop("'", deparse(substitute(module)), "' is not a 'network_module' object.")
  if(is.null(module$edges)) {
    warning("Argument 'module' contains no connections. Returning module unmodified.")
  } else {
    weights <- rdist(nrow(module$edges))
    module <- set_module_weights(module, weights)
  }

  return(module)
}



#' @inherit rewire_connections_to_node
#' @export
rewire_connections_to_node.network_module <- function(x,
                                                      node,
                                                      prob_rewire = 1,
                                                      weights = NULL,
                                                      alpha = 100,
                                                      beta = 1,
                                                      epsilon = 10^-5,
                                                      run_checks = TRUE,
                                                      ...) {
  if(!(class(x) == "network_module")) 
    stop("'", deparse(substitute(x)), "' is not a 'network_module' object.")
  
  adj_matrix <- get_adjacency_matrix(x)
  adj_matrix <- rewire_connections_to_node(adj_matrix, node, prob_rewire, 
                                           weights, alpha, beta, epsilon, run_checks, ...)
  x <- set_module_edges(x, adj_matrix)
  
  return(x)
}

#' @inherit rewire_connections
#' @export
rewire_connections.network_module <- function(x,
                                              prob_rewire = 1,
                                              weights = NULL,
                                              alpha = 100,
                                              beta = 1,
                                              epsilon = 10^-5,
                                              run_checks = TRUE, 
                                              ...) {
  nodes <- x$nodes
  x <- get_adjacency_matrix(x)
  x <- rewire_connections(x, prob_rewire, weights, alpha, beta, epsilon, run_checks, ...)
  x <- create_module_from_adjacency_matrix(x, nodes)
  return(x)
}

#' @inherit remove_connections_to_node
#' @export
remove_connections_to_node.network_module <- function(x,
                                                      node,
                                                      prob_remove,
                                                      run_checks = TRUE,
                                                      ...) {
  if(!(class(x) == "network_module")) 
    stop("'", deparse(substitute(x)), "' is not a 'network_module' object.")
  
  adj_matrix <- get_adjacency_matrix(x)
  adj_matrix <- remove_connections_to_node(adj_matrix, node, prob_remove, run_checks, ...)
  x <- set_module_edges(x, adj_matrix)
  
  return(x)
}

#' @inherit remove_connections
#' @export
remove_connections.network_module <- function(x,
                                              prob_remove,
                                              run_checks = TRUE,
                                              ...) {
  if(!(class(x) == "network_module")) 
    stop("'", deparse(substitute(x)), "' is not a 'network_module' object.")
  
  adj_matrix <- get_adjacency_matrix(x)
  adj_matrix <- remove_connections(adj_matrix, prob_remove, run_checks, ...)
  x <- set_module_edges(x, adj_matrix)
  
  return(x)
}



###########################################################################
#
# Utility functions for network_module.
#
###########################################################################

#' Print function for 'network_module' object.
#' 
#' @param x A 'network_module' object.
#' @param ... Additional arguments are ignored.
#' @return Prints a summary of the module.
#' @export
#' @examples 
#' module <- random_module(1:10)
#' module
#' print(module)
print.network_module <- function(x, ...) {
  if(!(class(x) == "network_module")) 
    stop("'", deparse(substitute(x)), "' is not a 'network_module' object.")
  
  vals <- get_network_characteristics(x)
  message <- paste0(ifelse(is_weighted(x), "A weighted", "An unweighted"), 
                    " module containing ", vals$p, " nodes and ", 
                    vals$n_edges, " edges.\n",
                    "Contains nodes: ",
                    ifelse(vals$p > 50, 
                           paste(paste(x$nodes[1:50], collapse = ", "), "..."),
                           paste(x$nodes, collapse = ", ")), "\n")
  
  cat(message)
  print(round(unlist(vals[-c(1, 2, 6)]), 3))
}

#' @inherit is_weighted
#' @export
is_weighted.network_module <- function(x, ...) {
  if(!(class(x) == "network_module")) 
    stop("'", deparse(substitute(x)), "' is not a 'network_module' object.")
  
  if(is.null(x$edges)) {
    return(TRUE)
  }
  if(ncol(x$edges) >= 3) {
    return(TRUE)
  } 
  return(FALSE)
}
