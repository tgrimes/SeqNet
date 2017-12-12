# Make class "network"
# Add: setClass(
#       class = "network",
#       slots = list(var_name = "numeric"),
#       validity = function(object) { #Check validity. }
#      )
# Add: setMethod(f = "initialize", signature = "network",
#                definition = function(.Object, arg1, arg2) {
#                 cat("--- bmiTraj: initializator --- \n")
#                 rownames(arg1) <- paste("I", 1:nrow(arg1), sep='')
#                 .Object@arg1 <- arg1
#                 .Object@arg2 <- arg2
#                 return(.Object)
#                }
#               )
# Add: setMethod("[", "network", function(object, i, j, drop) { })
# Add: setMethod("show", "network", function(object) { })
# Add: setMethod("print", "network", function(object) { })
# Add: setMethod("plot", "network", function(object) { })
# Add: setGeneric(name = "fun_name", def = function(object) standardGeneric("fun_name"))
#      setMethod(
#       f = "fun_name",
#       signature = "network",
#       definition = function(object) {}
#      )
#  with "fun" = "add_hub", "add_clique", etc.

setClass(Class = "network")

#' Create a network object.
#' 
#' Creates a graph with certain features. This network is then
#' used with other sample generating methods to obtain count data. 
#' @param p The number of nodes in the graph
#' @param n_cliques The number of cliques in the graph.
#' @param clique_size A vector containing the size of each clique.
#' @param n_hubs The number of hub genes in the graph.
#' @param hub_size A vector containing the size (hub degree + 1) of each hub.
#' @param cliques (optional) A list of vectors, with each vector specifying a
#'   set genes within a clique.
#' @param hubs (optional) A list of vectors, with each vector specifying a set
#'   of hub genes; the first entry in the vector is the hub.
#' @param nonoverlapping If true, cliques and hubs will not share any nodes.
#' @return A network object.
#' @export 
#' @examples 
#' create_network(p = 50, cliques = c(1:5, 6:10), hubs = c(11:15, 16:25))
create_network <- function(p, 
                           n_cliques = NULL, clique_size = NULL, cliques = NULL,
                           n_hubs = NULL, hub_size = NULL, hubs = NULL,
                           n_modules = NULL, module_size = NULL, modules = NULL,
                           nonoverlapping = FALSE, module_prob = 0.05) {

  if(p <= 0) stop("p should be positive.")
  
  # Updated when nonoverlapping is true to keep track of used nodes.
  node_set <- 1:p
  
  # Assign genes to structures.
  setup_structure <- function(struct_list, n_structs, struct_sizes) {
    # If all parameters are NULL, return NULL.
    if(is.null(struct_list) && is.null(n_structs) && is.null(struct_sizes)) {
      return(NULL)
    }
    
    # Set up structures.
    if(!is.null(struct_list)) {
      # Structures are provided.
      # Return `struct_list` as is.
    } else {
      # Structures are not provided; generate them based on `struct_sizes`.
      # If no struct_sizes are provided, generate random sizes based on n_structs.
      if(is.null(struct_sizes)) {
        struct_sizes <- 2 + round(rexp(n_structs, ifelse(p <= 1000, 10, 100) / p))
      }
      # If struct_sizes is of length 1, use this size for all n_structs structures.
      if(length(struct_sizes) == 1 && !is.null(n_structs)) {
        struct_sizes <- rep(struct_sizes, n_structs)
      }
      # If n_structs is not provided, set to length of struct_sizes.
      if(is.null(n_structs)) {
        n_structs <- length(struct_sizes)
      }
      # Check that length of struct_sizes is equal to n_structs.
      if(length(struct_sizes) != n_structs) {
        stop("Length of struct_size does not equal n_structs.")
      }
      
      struct_list <- vector("list", n_structs)
      if(n_structs > 0) {
        # Assign genes to each structure.
        for(i in 1:n_structs) {
          if((length(node_set) < struct_sizes[i]) && nonoverlapping) {
            stop("Number of nodes is too small for non-overlapping structures.")
          }
          struct_list[[i]] <- sample(node_set, struct_sizes[i], replace = FALSE)
          if(nonoverlapping) {
            node_set <<- setdiff(node_set, struct_list[[i]])
          }
        }
      }
    }
    
    return(struct_list)
  }
  
  # Set up cliques.
  cliques <- setup_structure(cliques, n_cliques, clique_size)
  cliques <- lapply(cliques, function(clique) {
    nodes <- clique
    struct <- rep(1, length(nodes))
    return(list(nodes = nodes,
                struct = struct))
  })
  
  # Set up hubs.
  hubs <- setup_structure(hubs, n_hubs, hub_size)
  hubs <- lapply(hubs, function(hub) {
    nodes <- hub
    struct <- rep(1, length(nodes))
    return(list(nodes = nodes,
                struct = struct))
  })
  
  # Set up modules.
  modules <- setup_structure(modules, n_modules, module_size)
  modules <- lapply(modules, function(module) {
    p <- length(module)
    nodes <- module
    # Create an adjacency matrix for the module.
    struct <- matrix(0, p, p) 
    # path <- sample(1:p) # Create a path through the module.
    # sapply(1:(p - 1), function(i) { struct[path[i], path[i + 1]] <<- 1 })
    
    # Add random structures to the module.
    struct <- struct + matrix(rbinom(length(nodes)^2, 1, module_prob), p, p)
    struct <- 1 * ((struct + t(struct)) > 0) # Symmetrize.
    diag(struct) <- 0
    return(list(nodes = nodes,
                struct = struct))
  })
  
  network <- list(p = p,
                  cliques = cliques, 
                  hubs = hubs, 
                  modules = modules, 
                  node_names = paste(1:p))
  
  # By default, add random sign to each connection in the network.
  network <- add_sign_to_network(network)
  
  class(network) <- "network"
  
  return(network)
}

#' Get adjacency matrix of a network
#' 
#' The adjacency matrix is constructed from all the structures in the network
#' and summarizes which nodes are connected. 
#' @param network The network to get adjacency matrix for.
#' @return A p by p adjacency matrix with entry ij = 1 if node i and j are 
#' connected, and 0 otherwise.
#' @export
get_adj_matrix_from_network <- function(network) {
  p <- network$p
  
  adj_matrix <- matrix(0, nrow = p, ncol = p)
  if(length(network$cliques) > 0) {
    lapply(network$cliques, function(clique) {
      adj_matrix[clique$nodes, clique$nodes] <<- 1
      diag(adj_matrix)[clique$nodes] <<- 0
    })
  }
  if(length(network$hubs) > 0) {
    lapply(network$hubs, function(hub) {
      adj_matrix[hub$nodes[1], hub$nodes[-1]] <<- 1
      adj_matrix[hub$nodes[-1], hub$nodes[1]] <<- 1
    })
  }
  if(length(network$modules) > 0) {
    lapply(network$modules, function(module) {
      nodes <- module$nodes
      adj_matrix[nodes, nodes] <<- 1 * (module$struct | adj_matrix[nodes, nodes])
    })
  }
  
  colnames(adj_matrix) <- network$node_names
  
  return(adj_matrix)
}
#' Add random signs to network edges
#' 
#' A random direction of association (positive or negative) is added to each 
#' connection in the cliques, hubs, and modules.
#' @param network The network to modify. Can be a network object or adjacency
#' matrix.
#' @return A network object with added structure to cliques, hubs, and modules.
#' @export
add_sign_to_network <- function(network) {
  network$hubs <- lapply(network$hubs, function(hub) {
    k <- length(hub$nodes) - 1
    hub$signs <- c(1, sample(c(-1, 1), k, replace = TRUE))
    return(hub)
  })
  
  network$cliques <- lapply(network$cliques, function(clique) {
    k <- length(clique$nodes)
    clique$signs <- sample(c(-1, 1), k, replace = TRUE)
    return(clique)
  })
  
  network$modules <- lapply(network$modules, function(module) {
    k <- length(module$nodes)
    module$signs <- matrix(sample(c(-1, 1), k^2, replace = TRUE), k, k)
    module$signs[upper.tri(module$signs)] <- 1
    module$signs <- module$signs * module$struct
    return(module)
  })
  
  return(network)
}


#' Add random weights to network edges
#' 
#' Used by the generator when generating samples; the strength of the connections
#' between genes varies from sample to sample. This function generates a random
#' strength for each edge from a truncated Normal distribution. The 
#' intensity parameter specifies the standard deviation of the Normal distribution.
#' @param network The network to modify.
#' @param intensity a positive number; this determines the scale of the weights.
#' @return A network with a weighted adjacency matrix. 
#' @export
add_weight_to_network <- function(network, intensity = 1) {
  adj_matrix <- network$adj_matrix
  p <- network$p
  weight_matrix <- matrix(0, p, p)
  
  lapply(network$hubs, function(hub) {
    nodes <- hub$nodes
    signs <- hub$signs[-1] # Determines pos. or neg. association with hub node.
    k <- length(nodes) - 1
    weights <- rnorm(k, 0, intensity) # Different weights for each edge.
    weight_matrix[nodes[1], nodes[-1]] <<- 
      weights * signs + weight_matrix[nodes[1], nodes[-1]]
    weight_matrix[nodes[-1], nodes[1]] <<- 
      weights + weight_matrix[nodes[-1], nodes[1]] 
  })
  
  lapply(network$cliques, function(clique) {
    nodes <- clique$nodes
    signs <- clique$signs # Determines pos. or neg. association with latent factor.
    k <- length(nodes)
    weight <- matrix(rnorm(1, 0, intensity / sqrt(k)), k, k) # Same weight for all edges.
    weight_matrix[nodes, nodes] <<- 
      weight %*% diag(signs) + weight_matrix[nodes, nodes]
  })
  
  lapply(network$modules, function(module) {
    nodes <- module$nodes
    struct <- module$struct
    signs <- module$struct
    k <- length(nodes)
    weight <- matrix(rnorm(k^2, 0, intensity / sqrt(2)), k, k) # New weight for each edge.
    weight <- weight + t(weight) # Symmetrize weights.
    weight <- weight * struct # Impose structure of the module.
    weight_matrix[nodes, nodes] <<- 
      weight * signs + weight_matrix[nodes, nodes]
  })
  
  network$weight_matrix <- weight_matrix
  
  return(network)
}


