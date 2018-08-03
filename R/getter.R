#' Get adjacency matrix of a network
#' 
#' The adjacency matrix is constructed from all the structures in the network
#' and summarizes which nodes are connected. 
#' @param network The network to get adjacency matrix for.
#' @param weighted If true, the ij'th entry will be the total number of 
#' connections between nodes i and j.
#' @return A p by p adjacency matrix with entry ij >= 1 if node i and j are 
#' connected, and 0 otherwise.
#' @export
get_adj_matrix_from_network <- function(network, weighted = FALSE) {
  p <- network$p
  
  A <- matrix(0, nrow = p, ncol = p) # Adjacency matrix.
  if(length(network$cliques) > 0) {
    lapply(network$cliques, function(clique) {
      A[clique$nodes, clique$nodes] <<- A[clique$nodes, clique$nodes] + 1
      diag(A)[clique$nodes] <<- 0
    })
  }
  if(length(network$hubs) > 0) {
    lapply(network$hubs, function(hub) {
      A[hub$nodes[1], hub$nodes[-1]] <<- A[hub$nodes[1], hub$nodes[-1]] + 1
      A[hub$nodes[-1], hub$nodes[1]] <<- A[hub$nodes[-1], hub$nodes[1]] + 1
    })
  }
  if(length(network$modules) > 0) {
    lapply(network$modules, function(module) {
      nodes <- module$nodes
      A[nodes, nodes] <<- A[nodes, nodes] + module$struct
    })
  }
  
  if(!weighted) {
    A <- 1 * (A > 0) # Reduce counts to 1 or 0.
  }
  
  colnames(A) <- network$node_names
  
  return(A)
}

#' Get summary for a node in the network.
#' 
#' @param network A network object.
#' @param node The node to summarize. Can be a character string if the nodes
#' are labeled, or a integer value from 1 to p.
#' @return A list containing summary, a list of strings detailing each 
#' structure the node belongs to; struct_count, the number of structures
#' the node belongs to; and degree, the total number of connections to the 
#' node (summed over all structures).
#' @export
get_summary_for_node <- function(network, node) {
  if(is.character(node)) {
    node <- which(network$node_names == node)
    if(length(node) == 0) {
      warning("Node not found amoung node names. Returning NULL.")
      return(NULL)
    }
  }
  
  degree <- 0
  struct_list <- vector("list", 0)
  struct_count <- 0
  add_struct <- function(string, struct_type = NULL) {
    struct_list[struct_count + 1] <<- string
    if(!is.null(struct_type)) {
      names(struct_list)[struct_count + 1] <<- struct_type
    }
    struct_count <<- struct_count + 1
  }
  
  if(length(network$hubs) > 0) {
    for(i in 1:length(network$hubs)) {
      hub_nodes <- network$hubs[[i]]$nodes
      node_index <- which(node == hub_nodes)
      if(length(node_index) == 1) {
        if(node_index[1] == 1) {
          n_connections <- length(hub_nodes)
          add_struct(paste0("Hub gene for hub ", i, ": ", 
                            n_connections, " connections"), 
                     "Hub")
          degree <- degree + n_connections
        } else {
          add_struct(paste0("Leaf gene in hub ", i), 
                     "Hub")
          degree <- degree + 1
        }
      }
    }
  }
  
  if(length(network$modules) > 0) {
    for(i in 1:length(network$modules)) {
      module_struct <- network$modules[[i]]$struct
      node_index <- which(node == network$modules[[i]]$nodes)
      if(length(node_index) == 1) {
        n_connections <- sum(module_struct[, node_index])
        add_struct(paste0("Gene in module ", i, ": ",
                          n_connections, " connections"), 
                   "Module")
        degree <- degree + n_connections
      }
    }
  }
  
  
  if(length(network$clique) > 0) {
    for(i in 1:length(network$clique)) {
      clique_nodes <- network$clique[[i]]$nodes
      node_index <- which(node == clique_nodes)
      if(length(node_index) == 1) {
        n_connections <- length(clique_nodes) - 1
        add_struct(paste0("Gene in clique ", i, ": ",
                          n_connections, " connections"), 
                   "Clique")
        degree <- degree + n_connections
      }
    }
  }
  
  return(list(summary = struct_list, 
              struct_count = struct_count, 
              degree = degree))
}

#' Get the degree distribution for a network.
#' 
#' Counts the connections to each node within each structure. Note, this
#' is not the same as the degree distribution from the adjacency matrix
#' obtained from the network, which collapses the individual structures into
#' one graph.
#' @param network A network object.
#' @return A vector of length p, containing the degree for each node in the 
#' network.
#' @export
get_degree_distribution <- function(network) {
  degree <- vector("numeric", network$p)
  for(i in 1:network$p) {
    degree[i] <- get_summary_for_node(network, i)$degree
  }
  return(degree)
}

#' Get the change in degree for each node between two networks.
#' 
#' Counts the connections to each node within each structure. Note, this
#' is not the same as the change in degree between the two adjacency matricies
#' of each network, which collapses individual structures.
#' @param network1 A network object for the first network.
#' @param network2 A network object for the second network.
#' @return A vector of length p, containing the degree change for each node in 
#' the two networks.
#' @export
get_degree_change <- function(network_1, network_2) {
  A1 <- get_adj_matrix_from_network(network_1, weighted = TRUE)
  A2 <- get_adj_matrix_from_network(network_2, weighted = TRUE)
  delta <- apply(A1 - A2, 2, sum) # Degree change.
  return(delta)
}

#' Get the clustering of nodes in the network
#' 
#' Structure is synonnymous with cluster or group. The structures in the 
#' network may overlap, meaning a node can be a member of multiple clusters. 
#' This function identifies which of the nodes belong to each structure.
#' @param network A network object.
#' @return A p by g matrix, M, where g is the number of structures in the network
#' and p is the number of nodes. M_ij = 1 indicates that node i is a member
#' of structure j. M_ij = 0 indicates non-membership.
#' @export
get_membership_matrix <- function(network) {
  n_hubs <- length(network$hubs)
  n_modules <- length(network$modules)
  n_cliques <- length(network$cliques)
  
  n_structures <- n_hubs + n_modules + n_cliques
  G <- matrix(0, nrow = network$p, ncol = n_structures)
  colnames(G) <- 1:n_structures
  
  if(n_hubs > 0) {
    for(i in 1:n_hubs) {
      G[network$hubs[[i]]$nodes, i] <- 1
      colnames(G)[i] <- paste0("hub_", i)
    }
  }
  
  if(n_modules > 0) {
    for(i in 1:n_modules) {
      if(any(apply(network$modules[[i]]$struct, 2, sum) == 0)) {
        warning("node in module but has no connections")
      }
      G[network$modules[[i]]$nodes, n_hubs + i] <- 1
      colnames(G)[n_hubs + i] <- paste0("module_", i)
    }
  }
  
  if(n_cliques > 0) {
    for(i in 1:n_cliques) {
      G[network$cliques[[i]]$nodes, n_hubs + n_modules + i] <- 1
      colnames(G)[n_hubs + n_modules + i] <- paste0("clique_", i)
    }
  }
  
  return(G)
}


#' Get a list of pathways in the network
#' 
#' Pathways refer to hubs, modules, and cliques. The pathways in the 
#' network may overlap, meaning a node can be a member of multiple pathways. 
#' This function identifies which of the nodes belong to each pathway.
#' @param network A network object.
#' @return A list of g vectors of integers, where g is the number of pathways 
#' in the network and each vector contains the indicies of the nodes contained
#' in the pathway.
#' @export
get_membership_list <- function(network) {
  n_hubs <- length(network$hubs)
  n_modules <- length(network$modules)
  n_cliques <- length(network$cliques)
  
  n_structures <- n_hubs + n_modules + n_cliques
  G <- vector("list", n_structures)

  if(n_hubs > 0) {
    for(i in 1:n_hubs) {
      G[[i]] <- network$hubs[[i]]$nodes
      names(G)[i] <- paste0("hub_", i)
    }
  }
  
  if(n_modules > 0) {
    for(i in 1:n_modules) {
      if(any(apply(network$modules[[i]]$struct, 2, sum) == 0)) {
        warning("node in module but has no connections")
      }
      G[[n_hubs + i]] <- network$modules[[i]]$nodes
      names(G)[n_hubs + i] <- paste0("module_", i)
    }
  }
  
  if(n_cliques > 0) {
    for(i in 1:n_cliques) {
      G[[n_hubs + n_modules + i]] <- network$cliques[[i]]$nodes
      names(G)[n_hubs + n_modules + i] <- paste0("clique_", i)
    }
  }
  
  return(G)
}