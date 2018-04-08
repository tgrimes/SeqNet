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
#' used with other sample generating methods to obtain count data. Note: the
#' module structure is used to incorporate general pathways in the graph. These
#' are randomly constructed by generating a small-world graph using the Watts-Strogatz
#' method (implemented through igraph::watts.strogatz.game). 
#' @param p The number of nodes in the graph
#' @param n_cliques The number of cliques in the graph.
#' @param clique_size A vector containing the size of each clique.
#' @param cliques A list of vectors, with each vector specifying a
#'   set genes within a clique.
#' @param n_hubs The number of hub genes in the graph.
#' @param hub_size A vector containing the size (hub degree + 1) of each hub.
#' @param hubs A list of vectors, with each vector specifying a set
#'   of hub genes; the first entry in the vector is the hub.
#' @param n_modules The number of modules in the graph.
#' @param module_size A vector containing the size of each hub.
#' @param modules A list of vectors, with each vector specifying a set of
#'   genes within each module.
#' @param nonoverlapping If true, cliques and hubs will not share any nodes.
#' @param module_neig Used for generating small-world networks for modules; 
#' specifies the initial degree of each node.
#' @param module_prob Used for generating small-world networks for modules; 
#' specifies the rewiring probability of each initial edges. At 0 the graph is
#' a lattice ring, and at 1 it is a random graph. 
#' @return A network object.
#' @export 
#' @examples 
#' # Create a network with no connections:
#' create_network(p = 100)
#' # Create a small-world network using all verticies:
#' create_network(p = 100, modules = list(1:100))
#' # Create a network with cliques and hubs:
#' create_network(p = 100, cliques = list(1:5, 6:10), hubs = list(11:15, 16:25))
create_network <- function(p, 
                           n_cliques = NULL, clique_size = NULL, cliques = NULL,
                           n_hubs = NULL, hub_size = NULL, hubs = NULL,
                           n_modules = NULL, module_size = NULL, modules = NULL,
                           nonoverlapping = FALSE, module_neig = 2, module_prob = 0.90) {

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
    struct <- rep(1L, length(nodes))
    return(list(nodes = nodes,
                struct = struct))
  })
  
  # Set up hubs.
  hubs <- setup_structure(hubs, n_hubs, hub_size)
  hubs <- lapply(hubs, function(hub) {
    nodes <- hub
    struct <- rep(1L, length(nodes))
    return(list(nodes = nodes,
                struct = struct))
  })
  
  # Set up modules.
  modules <- setup_structure(modules, n_modules, module_size)
  modules <- lapply(modules, function(module) {
    p <- length(module)
    nodes <- module
    # Generate a small world graph for the module using the Watts-Strogatz model.
    struct <- as.matrix(
                igraph::get.adjacency(
                  igraph::watts.strogatz.game(1, p, 
                                              module_neig, module_prob)))
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


#' Create a network object from adjacency matrix
#' 
#' @param adj_matrix The adjacency matrix for the network. This is converted
#' to a single module structure.
#' @param p The number of nodes in the network. Taken to be the number of nodes in
#' adjacency matrix by default.
#' @return A network object.
#' @export 
create_network_from_adjacency <- function(adj_matrix, p = ncol(adj_matrix)) {
  if(p < ncol(adj_matrix)) 
    stop("p must be >= ncol(adj_matrix).")
  if(!isSymmetric(unname(adj_matrix))) 
    stop("adj_matrix is not symmetric.")
  if(!all(adj_matrix %in% c(0, 1))) 
    stop("adj_matrix contains elements other than 0 or 1.")
  
  # Use node names from adj_matrix, if provided and no new nodes are added.
  if(!is.null(colnames(adj_matrix)) && p == ncol(adj_matrix)) {
    node_names <- colnames(adj_matrix)
  } else {
    node_names <- paste(1:p)
  }
  
  # Set up module.
  modules <- list(list(nodes = 1:ncol(adj_matrix),
                       struct = struct))
  network <- list(p = p,
                  cliques = NULL, 
                  hubs = NULL, 
                  modules = modules, 
                  node_names = node_names)
  
  # By default, add random sign to each connection in the network.
  network <- add_sign_to_network(network)
  
  class(network) <- "network"
  
  return(network)
}

create_network_scale_free <- function(p) {
  max_iter <- 500
  iter <- 1
  
  min_n_hubs <- max(5, floor(sqrt(p) / 4))
  min_n_modules <- max(1, floor(sqrt(p) / 4))
  max_n_hubs <- ceiling(sqrt(p) * 3)
  max_n_modules <- ceiling(sqrt(p) * 3)
  min_mu_hubs <- 2
  max_mu_hubs <- max(2, p * 0.25)
  min_mu_modules <- 2
  max_mu_modules <- max(2, p * 0.5)
  
  n_hubs <- rep(0, max_iter)
  n_modules <- rep(0, max_iter)
  mu_hubs <- rep(0, max_iter)
  mu_modules <- rep(0, max_iter)
  
  n_hubs[1] <- sample(min_n_hubs:max_n_hubs, 1) #rep(ceiling(sqrt(p) / 2), max_iter)
  n_modules[1] <- sample(min_n_modules:max_n_modules, 1) #rep(ceiling(sqrt(p) / 3), max_iter)
  mu_hubs[1] <- (min_mu_hubs + max_mu_hubs) / 2 #rep(sqrt(p), max_iter)
  mu_modules[1] <- (min_mu_modules + max_mu_modules) / 2 #rep(sqrt(p) * 4, max_iter)
  
  d_mu_hubs <- 0
  d_n_hubs <- 0
  d_mu_modules <- 0
  d_n_modules <- 0
  
  r_sq <- rep(0, max_iter)
  par(mfrow = c(1, 3))
  while(iter < max_iter) {
    network <- create_network(p = p,
                             hubs = lapply(rnbinom(ceiling(n_hubs[iter]), 2, 
                                                   mu = mu_hubs[iter]), function(x) {
                               x <- max(min(x, p), 2) # Make sure x is in [2, p].
                               sample(1:p, x)
                             }),
                             modules = lapply(rnbinom(ceiling(n_modules[iter]), 2, 
                                                      mu = mu_modules[iter]), function(x) {
                               x <- max(min(x, p), 2) # Make sure x is in [2, p].
                               sample(1:p, x)
                             }),
                             module_neig = 3,
                             module_prob = 0.9)
    
    # Compute r_sq for scale-free distribution of degree.
    degree <- apply(get_adj_matrix_from_network(network), 2, sum)
    #par(mfrow = c(1, 3))
    #r_sq <- WGCNA::scaleFreePlot(degree)$scaleFreeRsquared
    #hist(degree + 1)
    #hist(log(degree + 1))
    r_sq[iter] <- summary(lm(log(approx(density(degree + 1), xout = degree + 1)$y) ~ log(degree + 1)))$r.squared
    
    # Plot progress.
    
    
    # Walk in parameter space.
    # Make step size proportional to 1 - r_sq.
    step <- floor(10 * (1 - r_sq[iter]))
    # If r_sq improved, bias direction towards current direction.
    if((iter == 1) || (r_sq[iter - 1] <= r_sq[iter])) {
      direction <- 1
    } else {
      direction <- -1
    }
    
    d_n_hubs <- direction * ceiling(d_n_hubs * (1 - iter / max_iter)) + 
      sample(seq(-step, step), 1)
    d_n_modules <- direction * ceiling(d_n_modules * (1 - iter / max_iter)) + 
      sample(seq(-step, step), 1)
    d_mu_hubs <- direction * d_mu_hubs * (1 - iter / max_iter) + 
      rnorm(1, 0, mu_hubs[1] * 0.2) * (1 - r_sq[iter])
    d_mu_modules <- direction * d_mu_modules * (1 - iter / max_iter) +
      rnorm(1, 0, mu_modules[1] * 0.2) * (1 - r_sq[iter])
    # } else {
    #   d_n_hubs <- sample(seq(-step, step), 1)
    #   d_n_modules <- sample(seq(-step, step), 1)
    #   d_mu_hubs <- rnorm(1, 0, mu_hubs[1] * 0.2) * (1 - r_sq[iter])
    #   d_mu_modules <- rnorm(1, 0, mu_modules[1] * 0.2) * (1 - r_sq[iter])
    # }
    
    # cat("(", d_n_hubs, ", ", d_n_modules, ")\n", sep = "")
    
    # Update parameters.
    n_hubs[iter + 1] <-     min(max(min_n_hubs, n_hubs[iter] + d_n_hubs),
                                max_n_hubs)
    n_modules[iter + 1] <-  min(max(min_n_modules, n_modules[iter] + d_n_modules), 
                                max_n_modules)
    mu_hubs[iter + 1] <-    min(max(min_mu_hubs, mu_hubs[iter] + d_mu_hubs), 
                                max_mu_hubs)
    mu_modules[iter + 1] <- min(max(min_mu_modules, mu_modules[iter] + d_mu_modules), 
                                max_mu_modules)
    
    
    if(iter > 100 && r_sq[iter] > 0.9) {
      if((max(r_sq[(iter - 20):iter]) - min(r_sq[(iter - 20):iter])) < 0.1) {
        break
      }
    }
    
    iter <- iter + 1
  }
  plot(1:iter, r_sq[1:iter], type = "l", 
       xlab = "iter", xlim = c(1, max_iter), 
       ylab = "r_sq", ylim = c(0, 1))
  plot(n_hubs[1:iter], n_modules[1:iter], type = "l",
       xlab = "Number of hubs", xlim = c(min_n_hubs, max_n_hubs),
       ylab = "Number of modules", ylim = c(min_n_modules, max_n_modules))
  plot(mu_hubs[1:iter], mu_modules[1:iter], type = "l",
       xlab = "Average hub size", xlim = c(min_mu_hubs, max_mu_hubs),
       ylab = "Average modules size", ylim = c(min_mu_modules, max_mu_modules))
  mu_hubs <- mu_hubs[1:iter]
  mu_modules <- mu_modules[1:iter]
  n_hubs <- n_hubs[1:iter]
  n_modules <- n_modules[1:iter]
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
    hub$signs <- c(1L, sample(c(-1, 1), k, replace = TRUE))
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
    # weights <- rbeta(k, 5, 3) * intensity # New weight for each edge.
    weight_matrix[nodes[1], nodes[-1]] <<- 
      weights * signs + weight_matrix[nodes[1], nodes[-1]]
    weight_matrix[nodes[-1], nodes[1]] <<- 
      weights + weight_matrix[nodes[-1], nodes[1]] 
  })
  
  lapply(network$cliques, function(clique) {
    nodes <- clique$nodes
    signs <- clique$signs # Determines pos. or neg. association with latent factor.
    k <- length(nodes)
    weights <- matrix(rnorm(1, 0, intensity / sqrt(k)), k, k) # Same weight for all edges.
    # weights <- matrix(rbeta(1, 5, 3), k, k) * intensity / k # New weight for each edge.
    weight_matrix[nodes, nodes] <<- 
      weights %*% diag(signs) + weight_matrix[nodes, nodes]
  })
  
  lapply(network$modules, function(module) {
    nodes <- module$nodes
    struct <- module$struct
    signs <- module$struct
    k <- length(nodes)
    weights <- matrix(rnorm(k^2, 0, intensity / sqrt(2)), k, k) # New weight for each edge.
    # weights <- matrix(rbeta(k^2, 5, 3), k, k) * intensity / 2 # New weight for each edge.
    weights <- weights + t(weights) # Symmetrize weights.
    weights <- weights * struct # Impose structure of the module.
    weight_matrix[nodes, nodes] <<- 
      weights * signs + weight_matrix[nodes, nodes]
  })
  
  network$weight_matrix <- weight_matrix
  
  return(network)
}

#' Remove node from a hub in the network
#' 
#' @param network The network to modify.
#' @param hub_index The index of the hub in the network to modify.
#' @param node_index The index of the node in the hub to remove. If the index
#' is 1, which is the hub node, the entire hub is removed from the network.
#' @return The modified network.
#' @export
remove_node_from_hub <- function(network, hub_index, node_index) {
  if(node_index == 1) {
    # Node is hub gene itself; remove entire hub.
    network$hubs[hub_index] <- NULL
    return(network)
  }
  
  network$hubs[[hub_index]]$nodes <- 
    network$hubs[[hub_index]]$nodes[-node_index]
  network$hubs[[hub_index]]$struct <- 
    network$hubs[[hub_index]]$struct[-node_index]
  network$hubs[[hub_index]]$signs <- 
    network$hubs[[hub_index]]$signs[-node_index]
  return(network)
}

#' Remove node from a module in the network
#' 
#' @param network The network to modify.
#' @param module_index The index of the module in the network to modify.
#' @param node_index The index of the node in the hub to remove.
#' @return The modified network.
#' @export
remove_node_from_module <- function(network, module_index, node_index) {
  network$modules[[module_index]]$nodes <- 
    network$modules[[module_index]]$nodes[-node_index]
  network$modules[[module_index]]$struct <- 
    network$modules[[module_index]]$struct[-node_index, -node_index]
  network$modules[[module_index]]$signs <- 
    network$modules[[module_index]]$signs[-node_index, -node_index]
  return(network)
}

#' Remove node from a clique in the network
#' 
#' @param network The network to modify.
#' @param clique_index The index of the clique in the network to modify.
#' @param node_index The index of the node in the hub to remove.
#' @return The modified network.
#' @export
remove_node_from_clique <- function(network, clique_index, node_index) {
  network$cliques[[clique_index]]$nodes <- 
    network$cliques[[clique_index]]$nodes[-node_index]
  network$cliques[[clique_index]]$struct <- 
    network$cliques[[clique_index]]$struct[-node_index]
  network$cliques[[clique_index]]$signs <- 
    network$cliques[[clique_index]]$signs[-node_index]
  return(network)
}

#' Remove node from the network
#' 
#' @param network The network to modify.
#' @param node The node to remove. Can be a character string if the nodes
#' are labeled, or a integer value from 1 to p.
#' @return The modified network.
#' @export
remove_connections_to_node <- function(network, node) {
  if(is.character(node)) {
    node <- which(network$node_names == node)
    if(length(node) == 0) {
      warning("Node not found amoung node names. Returning network unchanged.")
      return(network)
    }
  }
  
  if(length(network$hubs) > 0) {
    hubs_removed <- 0
    for(i in 1:length(network$hubs)) {
      hub_index <- i - hubs_removed
      node_index <- which(node == network$hubs[[hub_index]]$nodes)
      if(length(node_index) == 1) {
        if(node_index == 1) hubs_removed <- hubs_removed + 1
        
        network <- remove_node_from_hub(network, hub_index, node_index)
      }
    }
  }
  
  if(length(network$modules) > 0) {
    for(i in 1:length(network$modules)) {
      node_index <- which(node == network$modules[[i]]$nodes)
      if(length(node_index) == 1) {
        network <- remove_node_from_module(network, i, node_index)
      }
    }
  }
  
  if(length(network$clique) > 0) {
    for(i in 1:length(network$clique)) {
      node_index <- which(node == network$clique[[i]]$nodes)
      if(length(node_index) == 1) {
        network <- remove_node_from_clique(network, i, node_index)
      }
    }
  }
  
  return(network)
}


#' Rewire node in the network.
#' 
#' Within each module the node is a member of, all connections are reset. 
#' New connections are made with a specified probablility.
#' @param network The network to modify.
#' @param node The node to rewire. Can be a character string if the nodes
#' are labeled, or a integer value from 1 to p.
#' @param p The probabiliy of an edge to be created.
#' @return The modified network.
#' @export
rewire_connections_to_node <- function(network, node, p = 0.2) {
  if(p > 1 || p < 0) stop("p should be between 0 and 1.")
  if(is.character(node)) {
    node <- which(network$node_names == node)
    if(length(node) == 0) {
      warning("Node not found amoung node names. Returning network unchanged.")
      return(network)
    }
  }
  
  
  if(length(network$modules) > 0) {
    for(i in 1:length(network$modules)) {
      node_index <- which(node == network$modules[[i]]$nodes)
      if(length(node_index) == 1) {
        m <- length(network$modules[[i]]$nodes)
        connections <- rbinom(m, 1, p)
        network$modules[[i]]$struct[node_index, ] <- connections
        network$modules[[i]]$struct[, node_index] <- connections
        network$modules[[i]]$signs[node_index, ] <- connections * (-1)^rbinom(m, 1, 0.5)
        network$modules[[i]]$signs[, node_index] <- connections * (-1)^rbinom(m, 1, 0.5)
      }
    }
  }
  
  return(network)
}