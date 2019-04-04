setClass(Class = "network")

#' Create a network object.
#' 
#' Creates a 'network' object containing no modules. 
#' @param p The number of nodes in the network
#' @return A network object.
#' @export 
create_empty_network <- function(p) {
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'p'.
  check_positive_integer(p, checklist)
  
  report_checks(checklist)
  ##################################
  
  network <- list(p = p,
                  modules = NULL, 
                  node_names = NULL)
  
  class(network) <- "network"
  
  return(network)
}


#' Create a network object.
#' 
#' Creates a graph with certain features. This network is then
#' used with other sample generating methods to obtain count data. Note: the
#' module structure is used to incorporate general pathways in the graph. These
#' are randomly constructed by generating a small-world graph using the Watts-Strogatz
#' method (implemented through igraph::watts.strogatz.game). 
#' @param p The number of nodes in the graph
#' @param modules A named list of network_module objects, or a list 
#' of numeric vectors. If numeric vectors are provided, modules
#' will be created using these vectors to specify the nodes in each module. 
#' @param node_names (optional) Vector of strings providing names for each node
#' in the graph. Default names are "1", "2", ..., "p".
#' @param ... Additional arguments passed to random_module(); only
#' used if modules need to be generated.
#' @return A network object.
#' @export 
create_network_from_modules <- function(p, 
                                        modules,
                                        node_names = as.character(1:p),
                                        ...) {
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'p'.
  check_positive_integer(p, checklist)
  
  # Check 'modules'.
  if(is.null(modules)) {
    ArgumentCheck::addError(
      msg = paste("Argument 'modules' must be a list of 'network_module'",
                  "objects or numeric vectors, or a single integer value",
                  "indicating the number of modules to generate."),
      argcheck = checklist
    )
  } else if(class(modules) == "list") {
    # Check each element in the list 'modules'.
    if(!all(sapply(modules, function(m) class(m) == "network_module")) &&
       !all(sapply(modules, is.numeric))) 
      ArgumentCheck::addError(
        msg = paste("Argument 'modules' must be a list of 'network_module'",
                    "objects or numeric vectors, or a single integer value",
                    "indicating the number of modules to generate."),
        argcheck = checklist
      )
  } else if(class(modules) == "network_module" || 
            (is.numeric(modules) && length(modules) > 1)) {
    # If 'modules' is provided but is not a list, correct without warning.
    modules <- list(modules)
  } else if(is.numeric(modules) && (modules %% 1 == 0)) {
    # If a single integer value is provided, generate random modules.
    modules <- populate_modules_for_network(modules, 
                                            p,
                                            ...)
  } else {
    ArgumentCheck::addError(
      msg = paste("Argument 'modules' must be a list of 'network_module'",
                  "objects or numeric vectors, or a single integer value",
                  "indicating the number of modules to generate."),
      argcheck = checklist
    )
  }
  
  report_checks(checklist)
  ##################################
  
  # Set up list of modules for the network.
  if(!is.null(modules) && all(sapply(modules, is.numeric))) {
    n_modules <- length(modules)
    module_list <- vector("list", n_modules)
    module_names <- names(modules)
    for(i in 1:n_modules) {
      module_list[[i]] <- random_module(nodes = modules[[i]],
                                        name = module_names[i],
                                        add_weights = FALSE,
                                        ...)
    }
    names(module_list) <- module_names
  } else {
    module_list <- modules
  }
  
  network <- create_empty_network(p)
  network <- add_modules_to_network(network, module_list)
  network <- set_node_names(network, node_names)
  
  return(network)
}


###########################################################################
#
# Functions to create network
#
###########################################################################

#' Create a network object from adjacency matrix
#' 
#' @param adj_matrix The adjacency matrix for the network. This is converted
#' to a single module structure.
#' @param ... Additional arguments passed to
#' create_module_from_adjacency_matrix(). 
#' @return A network object.
#' @export 
create_network_from_adjacency_matrix <- function(adjacency_matrix, ...) {
  p <- ncol(adjacency_matrix)
  
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'adjacency_matrix'
  check_adjacency_matrix(adjacency_matrix, checklist)
  
  report_checks(checklist)
  ##################################
  
  # Use node names from 'adjacency_matrix', if provided.
  if(!is.null(colnames(adjacency_matrix))) {
    node_names <- colnames(adjacency_matrix)
  } else {
    node_names <- paste(1:p)
  }
  
  # Set up module.
  module <- create_module_from_adjacency_matrix(adjacency_matrix,
                                                nodes = 1:p,
                                                ...)
  
  network <- create_network_from_modules(p, 
                                         modules = list(module), 
                                         node_names = node_names)
  
  return(network)
}

#' Create a network object from an association matrix
#' 
#' @param association_matrix The association matrix for the network. This is converted
#' to a single module structure with partial correlations specified by the
#' matrix.
#' @param standardize If True, then 'association_matrix' will be standardized
#' to obtain partial correlations.
#' @param ... Additional arguments passed to 
#' create_module_from_association_matrix().
#' @return A network object.
#' @export 
create_network_from_association_matrix <- function(association_matrix, 
                                                   standardize = TRUE,
                                                   ...) {
  p <- ncol(association_matrix)
  
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'association_matrix'
  check_precision_matrix(association_matrix, checklist)
  
  report_checks(checklist)
  ##################################
  
  if(!all(diag(association_matrix) == 1) && standardize) {
    cat("Standardizing association_matrix to obtain partial correlations.\n")
    association_matrix <- cov2cor(association_matrix)
  }
  
  # Use node names from association_matrix, if provided and no new nodes are added.
  if(!is.null(colnames(association_matrix))) {
    node_names <- colnames(association_matrix)
  } else {
    node_names <- paste(1:p)
  }
  
  # Set up module.
  module <- create_module_from_association_matrix(association_matrix,
                                                  nodes = 1:p,
                                                  ...)
  
  network <- create_network_from_modules(p, 
                                         modules = list(module), 
                                         node_names = node_names)
  
  return(network)
}


#' Create a network object.
#' 
#' Creates an unweighted 'network' object containing randomly generated 
#' modules.
#' @param p The number of nodes in the network; 'p' is required to be 
#' between 10 and 20000.
#' @param n_modules The number of modules to include in the network. The
#' default value is 5\% of the number of nodes in the network.
#' @param min_module_size The minimum number of nodes in a module.
#' @param max_module_size (optional) If specified, any generated module sizes 
#' above this value will be reduced to 'max_module_size'.
#' @param avg_module_size The average number of nodes in a module.
#' @param sd_module_size The standard deviation of module size.
#' @param repeated_module_factor A value between 0 and 1. A gene selected for a 
#' module has its probability of being selected for another module multiplied by 
#' this factor.
#' @param ... Additional arguments passed to 'random_module()' and, further 
#' downstream, to 'update_module_with_random_struct()'; those of particular 
#' interest may include the 'lattice_neig' and 'rewire_prob' arguments.
#' @return An unweighted network object.
#' @export 
random_network <- function(p,
                           n_modules = ceiling(p * 0.05),
                           min_module_size = 10,
                           max_module_size = Inf,
                           avg_module_size = 25,
                           sd_module_size = 20,
                           repeated_module_factor = 0.2,
                           ...) {
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'p'.
  check_positive_integer(p, checklist)
  check_in_closed_interval(p, checklist, 10, 20000)
  
  # Check 'n_modules'.
  check_positive_integer(n_modules, checklist)
  
  report_checks(checklist)
  ##################################
  
  nodes_list <- populate_modules_for_network(n_modules, 
                                             p,
                                             min_module_size,
                                             max_module_size,
                                             avg_module_size,
                                             sd_module_size, 
                                             repeated_module_factor)
  
  # Create the modules. 
  module_list <- lapply(nodes_list, random_module, add_weights = FALSE, ...)
  network <- create_network_from_modules(p, modules = module_list)
  network <- trim_modules(network)
  
  return(network)
}


#' Randomly sample subsets of genes for each module
#' 
#' Creates a collection of modules containing randomly samples genes.
#' @param n_modules The number of modules to include in the network.
#' @param p The number of nodes in the network.
#' @param min_module_size The minimum number of nodes in a module.
#' @param max_module_size (optional) If specified, any generated module sizes 
#' above this value will be reduced to 'max_module_size'.
#' @param avg_module_size The average number of nodes in a module.
#' @param sd_module_size The standard deviation of module size.
#' @param repeated_module_factor A value between 0 and 1. A gene selected for a 
#' module has its probability of being selected for another module multiplied by 
#' this factor.
#' @return A list containing the indicies for genes contained in each module.
#' @export 
populate_modules_for_network <- function(n_modules, 
                                         p,
                                         min_module_size = 10,
                                         max_module_size = Inf,
                                         avg_module_size = 25,
                                         sd_module_size = 10,
                                         repeated_module_factor = 0.2,
                                         ...) {
  if(repeated_module_factor < 0 || repeated_module_factor > 1)
    stop("Argument 'repeated_module_factor' must be between 0 and 1.")

  mu = (avg_module_size - min_module_size)
  size <- mu^2 / (sd_module_size^2 - mu)
  
  # Generate random size for each module.
  module_sizes <- min_module_size + 
    rnbinom(n_modules, size = size, mu = mu)
  # Reduce any sizes that are above 'max_module_size'.
  sizes_above_max <- which(module_sizes > max_module_size)
  if(length(sizes_above_max) > 0) {
    module_sizes[sizes_above_max] <- max_module_size
  }
  
  # Sample nodes to populate each module.
  all_nodes <- 1:p
  prob <- rep(1 / p, p)
  nodes_list <- vector("list", n_modules) 
  nodes_selected <- NULL
  for(i in 1:n_modules) {
    m <- module_sizes[i]
    nodes <- sort(sample(all_nodes, min(m, p), prob = prob))
    nodes_selected <- union(nodes, nodes_selected)
    prob <- rep(1 / p, p)
    prob[nodes_selected] <- prob[nodes_selected] * repeated_module_factor
    prob <- prob / sum(prob)
    nodes_list[[i]] <- nodes
  }
  
  return(nodes_list)
}

###########################################################################
#
# Getter functions for network objects.
#
###########################################################################

#' Get adjacency matrix of a network
#' 
#' The adjacency matrix is constructed from all modules in the network
#' and summarizes the connectivity among nodes. 
#' @param network A 'network' object; can be either weighted or unweighted.
#' @return An adjacency matrix with entry ij = 1 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_adjacency_matrix.network <- function(network, ...) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  if(is_weighted(network)) {
    adj_matrix <- get_association_matrix(network)
    adj_matrix <- 1 * (adj_matrix != 0)
  } else {
    p <- network$p
    adj_matrix <- matrix(0, nrow = p, ncol = p) # Adjacency matrix.
    n_modules <- length(network$modules)
    if(n_modules > 0) {
      for(i in 1:n_modules) {
        module_index <- network$modules[[i]]$nodes
        module_matrix <- get_adjacency_matrix.network_module(network$modules[[i]])
        adj_matrix[module_index, module_index] <- 
          adj_matrix[module_index, module_index] + module_matrix
      }
    }
    
    adj_matrix <- 1 * (adj_matrix > 0) # Reduce counts to 1 or 0.
  }
  
  
  colnames(adj_matrix) <- network$node_names
  
  return(adj_matrix)
}
# setMethod(
#   f = "get_adjacency_matrix",
#   signature = "network",
#   definition = get_adjacency_matrix.network
# )

#' Get the association matrix of a network
#' 
#' The association matrix for a weighted network is recovered from its
#' covariance structure. The association matrix will contain ones along the
#' diagonal, and off-diagonal entries are interpreted as partial correlations.
#' @param network A weighted 'network' object.
#' @return An association matrix with entry ij != 0 if node i and j are 
#' connected, and 0 otherwise. 
#' @export
get_association_matrix.network <- function(network, ...) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  if(!all(sapply(network$modules, is_weighted))) 
    stop("Network modules must all have weighted connections.")
  
  if(network$p > 5000)
    stop("The network is too large.")
  
  sigma <- get_sigma(network)
  assoc_matrix <- -solve(sigma)
  diag(assoc_matrix) <- -diag(assoc_matrix)
  assoc_matrix <- cov2cor(assoc_matrix)
  
  # Set values near zero to exactly zero.
  assoc_matrix[abs(assoc_matrix) <= 10^-13] <- 0
  diag(assoc_matrix) <- 0
  
  colnames(assoc_matrix) <- network$node_names
  
  return(assoc_matrix)
}

#' Get the covariance matrix of a network
#' 
#' The associations in each module are taken as partial correlations, and
#' the covariance matrix is calculated from these assuming that expression
#' for gene i is the weighted average over each module using 1/sqrt(m_i) 
#' as the weight, where m_i is the number of modules containing gene i.
#' @param network A weighted 'network' object.
#' @return A covariance matrix for the network.
#' @export
get_sigma.network <- function(network, ...) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  if(!all(sapply(network$modules, is_weighted))) 
    stop("Network modules must all have weighted connections.")
  
  m <- rep(0, network$p)
  sigma <- diag(0, network$p)
  for(i in 1:length(network$modules)) {
    nodes <- network$modules[[i]]$nodes
    sigma[nodes, nodes] <- sigma[nodes, nodes] + get_sigma(network$modules[[i]])
    m[nodes] <- m[nodes] + 1
  }
  
  index <- which(m == 0)
  m[index] <- 1
  diag(sigma)[index] <- 1
  # Cov of x_i and x_j is divided by by 1 / sqrt(m_i * m_j), where
  # m_i and m_j are the number of modules genes i and j belong to, respectively.
  # This assumes the data are generated by summing the expression of each
  # module and dividing by 1 / sqrt(m).
  denom <- 1 / sqrt(m)
  sigma <- diag(denom) %*% sigma %*% diag(denom)
  
  return(sigma)
}

#' Get summary for a node in the network.
#' 
#' @param node The node to summarize. Can be a character string 
#' corresponding to a name of a node in the network, or an integer value from 
#' 1 to p corresponding to the index of a node.
#' @param network A network object.
#' @return A list containing summary, a list of strings detailing each 
#' structure the node belongs to; struct_count, the number of structures
#' the node belongs to; and degree, the total number of connections to the 
#' node (summed over all structures).
#' @export
get_summary_for_node <- function(node, network) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'node'.
  if(is.character(node)) {
    # If a string is provided, find the index for this node in the network.
    node_index <- which(network$node_names == node)
    if(length(node) == 0) 
      ArgumentCheck::addError(
        msg = paste0("Argument 'network' does not contain a node nammed ", '"', node, '".'),
        argcheck = checklist
      )
  } else if(is.numeric(node) && (node %% 1 == 0)) {
    if(node < 1 || node > network$p) 
      ArgumentCheck::addError(
        msg = paste0("Argument 'node' = ", node, " must index a node in the network;",
                     " this network contains 1:", network$p, " nodes."),
        argcheck = checklist
      )
    node_index <- node
  } else {
    ArgumentCheck::addError(
      msg = "Argument 'node' must be a character string or integer value.",
      argcheck = checklist
    )
  }
  
  report_checks(checklist)
  ##################################
  
  # Initialize summary variables for node.
  degree <- 0
  n_modules <- 0
  n_connections <- 0
  module_index <- NULL
  
  # Update degree using the marginal adjacency matrix for the network.
  adj_matrix <- get_adjacency_matrix(network)
  degree <- sum(adj_matrix[, node_index])
  
  # Update module information.
  if(length(network$modules) > 0) {
    # Loop through each module.
    for(i in 1:length(network$modules)) {
      # If node is in this module, update summary information.
      if(node_index %in% network$modules[[i]]$nodes) {
        n_modules <- n_modules + 1
        module_index <- c(module_index, i)
        node_index_in_module <- which(network$modules[[i]]$nodes == node_index)
        adj_matrix <- get_adjacency_matrix(network$modules[[i]])
        n_connections <- n_connections + sum(adj_matrix[, node_index_in_module])
      }
    }
  }
  
  
  return(list(node = node,
              node_index = node_index,
              degree = degree,
              n_modules = n_modules,
              n_connections = n_connections,
              module_index = module_index))
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
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  adj_matrix <- get_adjacency_matrix(network)
  degree <- apply(adj_matrix, 2, sum)
  return(degree)
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
get_pathway_list <- function(network) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
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


#' Get node names of a network
#' 
#' @param network The 'network' object to get node names from.
#' @return A vector containing the node names
#' @export
get_node_names.network <- function(network, ...) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  return(network$node_names)
}

#' Get a list of modules from the network
#' 
#' @param network A 'network' object.
#' @return A list whose length is the number of modules in the network; 
#' each element is a vector containing the indicies of the nodes
#' that belong to that module.
#' @export
get_network_modules <- function(network) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  lapply(network$modules, function(module) module$nodes)
}

###########################################################################
#
# Utility functions for network objects.
#
###########################################################################

#' Adds a set of modules to the network
#' 
#' @param network The network to modify.
#' @param module_list A list of 'network_module' objects to add to the network.
#' @return The modified network.
#' @export
add_modules_to_network <- function(network, module_list) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'module_list'
  if(class(module_list) == "network_module") {
    # A single module is provided; put into a list without warning.
    module_list <- list(module_list)
  }
  # Only perform checks if module_list is not NULL and non-empty.
  if(!is.null(module_list) && length(module_list) > 0) {
    if(!all(sapply(module_list, function(m) class(m) == "network_module"))) { 
      ArgumentCheck::addError(
        msg = "Argument 'module_list' must be a list of 'network_module' objects.",
        argcheck = checklist
      )
    } else {
      # All elements in 'module_list' are modules; check that their nodes 
      # index nodes within the network.
      network_nodes <- 1:network$p
      index_bad_modules <- which(sapply(module_list, function(module) {
        !all(module$nodes %in% network_nodes)
      }))
      if(length(index_bad_modules) > 0)
        ArgumentCheck::addError(
          msg = paste0("'module_list' contains modules with nodes outside of",
                       " 1:", network$p, ". Errors with modules: ", 
                       paste(index_bad_modules, collapse = ", ")),
          argcheck = checklist
        )
    }
  }
  
  
  report_checks(checklist)
  ##################################
  
  module_names <- names(module_list)
  if(!is.null(module_names)) {
    for(i in 1:length(module_list)) {
      module_list[[i]] <- set_module_name(module_list[[i]], module_names[i])
    }
  }
  network$modules <- c(network$modules, module_list)
  
  return(network)
}


#' Replaces a module in the network
#' 
#' @param module_index The index of the module to replace.
#' @param module The new module to replace with.
#' @param network The network to modify.
#' @return The modified network.
#' @export
replace_module_in_network <- function(module_index, module, network) {
  if(!(class(module) == "network_module")) 
    stop(paste0("'", deparse(substitute(module)), 
                "' is not a 'network_module' object."))
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'module_index'.
  check_positive_integer(module_index, checklist)
  if(!(module_index %in% 1:length(network$modules)))
    ArgumentCheck::addError(
      msg = paste0("Argument 'module_index' must be in 1:", length(network$modules), "."),
      argcheck = checklist
    )
  
  # Check 'module'.
  if(!all(module$nodes %in% 1:network$p))
    ArgumentCheck::addError(
      msg = paste0("Argument 'module' must contain nodes in 1:", network$p, "."),
      argcheck = checklist
    )
  
  report_checks(checklist)
  ##################################
  
  network$modules[[module_index]] <- module
  
  return(network)
}


#' Removes the weights of all connections in the network
#' 
#' @param network The 'network' object to modify.
#' @return The modified 'network' object.
#' @export
remove_weights.network <- function(network) {
  network$modules <- lapply(network$modules, remove_weights)
  return(network)
}

#' Adds a random module of a given size to the network
#' 
#' @param network The 'network' object to modify.
#' @param module_size The size of the module to generate.
#' @param ... Additional arguments passed into random_module().
#' @return The modified 'network' object.
#' @export
add_random_module_to_network <- function(network, module_size, ...) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'module_size'
  if(module_size > network$p) 
    ArgumentCheck::addError(
      msg = paste0("Argument 'module_size' = ", module_size, " cannot be greater than the",
                   " number of nodes in the network (p =", network$p, ")"),
      argcheck = checklist
    )
  
  report_checks(checklist)
  ##################################
  
  module <- random_module(sample(1:network$p, module_size), ...)
  network <- add_modules_to_network(network, module)
  
  return(network)
}


#' Set the node names in a network
#' 
#' @param network The network to modify.
#' @param node_names A vector of strings containing the names for each node
#' in the network. If a numeric vector is provided, the values will be coerced
#' into strings. If 'node_names' is NULL, then the names will default to 
#' "1", "2", ..., "p".
#' @return The modified network.
#' @export
set_node_names <- function(network, node_names) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'node_names'
  if(is.null(node_names)) {
    # Use default names.
    node_names <- as.character(1:network$p)
    ArgumentCheck::addMessage(
      msg = paste0("Argument 'node_names' is NULL; setting default names to ",
                   '"1", "2", ..., "', network$p, '".'),
      argcheck = checklist
    )
  }
  if(is.numeric(node_names)) {
    # Coerce numeric values to character strings without warning.
    node_names <- as.character(node_names)
  }
  if(!is.character(node_names)) 
    ArgumentCheck::addError(
      msg = "Argument 'node_names' must be a vector of strings (or numeric vector).",
      argcheck = checklist
    )
  if(length(node_names) != network$p)
    ArgumentCheck::addError(
      msg = paste0("Length of argument 'node_names' must equal the number of nodes",
                   " in 'network' (p = ", network$p, ")."),
      argcheck = checklist
    )
  if(length(unique(node_names)) != length(node_names))
    ArgumentCheck::addError(
      msg = paste0("Argument 'node_names' must contain unique names."),
      argcheck = checklist
    )
    
  report_checks(checklist)
  ##################################
  
  network$node_names <- node_names
  
  return(network)
}

#' Rewire node in the network.
#' 
#' Within each module the node is a member of, all connections are reset. 
#' New connections are made with a specified probablility.
#' @param node The node to rewire. Can be a character string 
#' corresponding to a name of a node in the network, or an integer value from 
#' 1 to p corresponding to the index of a node.
#' @param network The network to modify.
#' @param rewire_prob The node at 'node_index' has probaility equal to 
#' 'rewire_prob' to be connected with any other node in the module.
#' @return The modified network.
#' @export
rewire_connections_to_node <- function(node, network, rewire_prob = 0.2) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  # If a vector of nodes are provided, loop over each one.
  if(length(node) > 1) {
    if(length(rewire_prob) == 1) {
      rewire_prob <- rep(rewire_prob, length(node))
    } else if(length(rewire_prob) != length(node)) {
      stop("Length of argument 'node' and 'rewire_prob' must be of the same length.")
    }
    for(i in 1:length(node)) {
      network <- rewire_connections_to_node(node[i], network, rewire_prob[i])
    }
    return(network)
  }
  
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'node'.
  if(is.character(node)) {
    # If a string is provided, find the index for this node in the network.
    node <- which(network$node_names == node)
    if(length(node) == 0) 
      ArgumentCheck::addError(
        msg = paste0("Argument 'network' does not contain a node nammed ", '"', node, '".'),
        argcheck = checklist
      )
  } else if(is.numeric(node) && (node %% 1 == 0)) {
    if(node < 1 || node > network$p) 
      ArgumentCheck::addError(
        msg = paste0("Argument 'node' = ", node, " must index a node in the network;",
                     " this network contains 1:", network$p, " nodes."),
        argcheck = checklist
      )
  } else {
    ArgumentCheck::addError(
      msg = "Argument 'node' must be a character string or integer value.",
      argcheck = checklist
    )
  }
  
  # Check 'rewire_prob'. (Allowed to equal 0 or 1.)
  check_in_closed_interval(rewire_prob, checklist, 0, 1)
  
  report_checks(checklist)
  ##################################
  
  if(length(network$modules) > 0) {
    for(i in 1:length(network$modules)) {
      if(node %in%  network$modules[[i]]$nodes) {
        network$modules[[i]] <- 
          rewire_connections_to_node_in_module(node,
                                               network$modules[[i]], 
                                               rewire_prob)
      }
    }
  } else {
    warning("Argument 'network' contains no modules. Returning network unmodified.")
  }
  
  return(network)
}


#' Remove all connections to a node
#' 
#' @param node The node to cut connections from. Can be a character string 
#' corresponding to a name of a node in the network, or an integer value from 
#' 1 to p corresponding to the index of a node.
#' @param network The network to modify.
#' @return The modified network.
#' @export
remove_connections_to_node <- function(node, network) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  return(rewire_connections_to_node(node, network, rewire_prob = 0))
}


#' Subsets all network modules to their largest components
#' 
#' Any nodes that are disconnected from the largest component in each module
#' are removed from that module. No nodes are removed from the network itself.
#' @param network The network to modify.
#' @return The modified network.
#' @export
trim_modules <- function(network) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  n_modules <- length(network$modules)
  if(n_modules == 0) {
    warning("Argument 'network' contains no modules. Returning network unmodified.")
  } else {
    for(i in 1:n_modules) {
      module <- network$modules[[i]]
      new_module <- remove_small_components_from_module(module)
      network <- replace_module_in_network(i, new_module, network)
    }
  }
  
  return(network)
}


#' Perturbs the connections in a network
#' 
#' By default, hub nodes have a 50% chance of being turned "off", and all other
#' nodes have a 10% chance of being perturbed by rewiring. No new hub nodes are
#' be created.
#' @param network The network to modify.
#' @param rewire_hub_prob The probability that any given hub node will have its
#' connections removed.
#' @param rewire_other_prob The probability that any given non-hub node will be
#' perturbed by rewiring. If a node is rewired, it will have, on average, the
#' same number of connections as its original state.
#' @return The modified network.
#' @export
perturb_network <- function(network, 
                            rewire_hub_prob = 0.5,
                            rewire_other_prob = 0.1) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  # Hub genes will be identified by having a degree above 3 SDs from the mean.
  degrees <- get_degree_distribution(network)
  index_hubs <- which(degrees > floor(mean(degrees) + 3 * sd(degrees)))
  index_others <- setdiff(1:network$p, index_hubs)
  rewire_prob <- degrees / network$p
  # Rewire hub nodes.
  if(length(index_hubs) > 0) {
    rewire <- rbinom(length(index_hubs), 1, rewire_hub_prob) == 1
    for(i in index_hubs[rewire]) {
      network <- remove_connections_to_node(i, network)
    }
  }
  # Rewire all other nodes.
  if(length(index_others) > 0) {
    rewire <- rbinom(length(index_others), 1, rewire_other_prob) == 1
    for(i in index_others[rewire]) {
      network <- rewire_connections_to_node(i, network, rewire_prob[i])
    }
  }
  
  return(network)
}


#' Print function for 'network' object.
#' 
#' @param network A 'network' object.
#' @param ... Additional arguments are ignored.
#' @return Prints a summary of the module.
#' @export
print.network <- function(network, ...) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  p <- network$p
  m <- length(network$modules)
  message <- paste0("Network containing ", p, " nodes and ", m, 
                    " module", ifelse(m > 1, "s", ""), ".\n")
  cat(message)
}


#' Check if a network is weighted
#' 
#' @param network The 'network' object to check.
#' @return A boolean value that is TRUE if all of the modules in the network
#' are weighted and FALSE otherwise. If the network contains no modules or
#' no connections, then this function returns TRUE.
#' @export
is_weighted.network <- function(network, ...) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  if(length(network$modules) == 0) {
    return(TRUE)
  } 
  
  if(all(sapply(network$modules, is_weighted))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Collapses all modules in network into a single module
#' 
#' This modification can be used if it is desired to simulate from a single
#' GGM rather than averaging over the GGMs for each module.
#' @param network The 'network' object to modify
#' @return The modified 'network' object.
#' @export
as_single_module <- function(network) {
  if(is_weighted(network)) {
    warning("A weighted network was provided; all connections weights have been removed.")
    network <- remove_weights(network)
  }
  network <- create_network_from_adjacency_matrix(
    get_adjacency_matrix(network))
  
  return(network)
}