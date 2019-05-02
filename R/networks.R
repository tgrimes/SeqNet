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
#' @param module_list A named list of 'network_module' objects.
#' @param node_names (optional) Vector of strings providing names for each node
#' in the graph. Default names are "1", "2", ..., "p".
#' @param ... Additional arguments passed to random_module(); only
#' used if modules need to be generated.
#' @return A network object.
#' @export 
create_network_from_modules <- function(p, 
                                        module_list, 
                                        node_names = as.character(1:p), 
                                        ...) {
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'p'.
  check_positive_integer(p, checklist)
  
  # Check 'module_list'.
  if(is.null(module_list)) {
    # NULL module_list is ok.
  } else if(class(module_list) == "list") {
    # Check each element in the list 'modules'.
    if(!all(sapply(module_list, function(m) class(m) == "network_module"))) 
      ArgumentCheck::addError(
        msg = paste("Argument 'module_list' must be a list of 'network_module'."),
        argcheck = checklist
      )
  } else if(class(module_list) == "network_module") {
    # If 'module_list' is provided but is not a list, correct without warning.
    module_list <- list(module_list)
  } else {
    ArgumentCheck::addError(
      msg = paste("Argument 'module_list' must be a list of 'network_module'."),
      argcheck = checklist
    )
  }
  
  report_checks(checklist)
  ##################################
  
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
                                         module_list = list(module), 
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
  if(all(diag(association_matrix) == 0)) {
    diag(association_matrix) <- 1
  }
  
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
                                         module_list = list(module), 
                                         node_names = node_names)
  
  return(network)
}


#' Create a network object.
#' 
#' Creates an unweighted 'network' object containing randomly generated 
#' modules.
#' @param p The number of nodes in the network; 'p' is required to be 
#' between 10 and 20000.
#' @param n_modules The number of modules to include in the network. If NULL,
#' then modules are created until all nodes in the network have positive degree.
#' @param consistent_connections If TRUE, then each module is modified so that,
#' if two genes are connected in one module, then they are connected in 
#' every module.
#' @param ... Additional arguments passed to 'create_modules_for_network()'.
#' @return An unweighted network object.
#' @export 
random_network <- function(p,
                           n_modules = NULL,
                           consistent_connections = TRUE,
                           ...) {
  ##################################
  # Check arguments for errors.
  checklist <- new_checklist()
  
  # Check 'p'.
  check_positive_integer(p, checklist)
  check_in_closed_interval(p, checklist, 10, 20000)
  
  # Check 'n_modules', if it is not NULL.
  check_nonnegative_integer(p, checklist)
  
  report_checks(checklist)
  ##################################
  
  module_list <- create_modules_for_network(n_modules, p, ...)
  network <- create_network_from_modules(p, module_list = module_list)
  
  # If two genes are connected in one module, make them connected in all modules.
  if(consistent_connections) {
    adj <- get_adjacency_matrix(network)
    module_list = lapply(network$modules, function(m) {
                           create_module_from_adjacency_matrix(adj[m$nodes, m$nodes], 
                                                               m$nodes)
                         })
    network <- 
      create_network_from_modules(p, module_list)
  }
  
  return(network)
}


#' Randomly sample subsets of genes for each module
#' 
#' Creates a collection of modules containing randomly samples genes.
#' @param n_modules The number of modules to include in the network.
#' @param p The number of nodes in the network.
#' @param avg_module_size The average number of nodes in a module.
#' @param sd_module_size The standard deviation of module size.
#' @param min_module_size The minimum number of nodes in a module.
#' @param max_module_size A positive value. Any generated module sizes above this 
#' value will be reduced to 'max_module_size'. Set to 'Inf' to avoid this 
#' truncation.
#' @param selection_weight A positive value used for sampling nodes for a new 
#' module.
#' @param ... Additional arguments passed to random_module().
#' @return A list containing the indicies for genes contained in each module.
#' @export 
create_modules_for_network <- function(n_modules, 
                                       p, 
                                       avg_module_size = 12, 
                                       sd_module_size = 4,
                                       min_module_size = 10, 
                                       max_module_size = 15, 
                                       selection_weight = 100,
                                       ...) {
  ############################################################
  # Check parameters for generating module size.
  ############################################################
  if(n_modules == 0) {
    return(NULL)
  }
  
  if(is.null(avg_module_size)) {
    avg_module_size <- 15
  }
  if(is.null(sd_module_size)) {
    sd_module_size <- avg_module_size / 2
  }
  
  if(avg_module_size <= min_module_size) {
    avg_module_size <- 1.1 * min_module_size
    warning(paste("Argument min_module_size should be smaller than avg_module_size.",
                  "Setting avg_module_size to 1.1 * min_module_size =", avg_module_size))
  }
  
  mu = (avg_module_size - min_module_size)
  if((sd_module_size^2 - mu) <= 0) {
    sd_module_size <- sqrt(mu) * 1.01
    warning(paste("Argument 'sd_module_size' is too small.",
                  "Setting to 1.01 * (avg_module_size - min_module_size) =", sd_module_size))
  }
  size <- mu^2/(sd_module_size^2 - mu)
  max_module_size <- min(p, max_module_size)
  
  ############################################################
  # Generate module sizes.
  ############################################################
  if(is.null(n_modules)) {
    # Generate extra modules sizes. Should not need more than p.
    module_sizes <- min_module_size + rnbinom(p, size = size, 
                                              mu = mu)
  } else {
    module_sizes <- min_module_size + rnbinom(n_modules, size = size, 
                                              mu = mu)
  }
  module_sizes <- sapply(module_sizes, min, max_module_size)
  
  ############################################################
  # Initialize nodes in the network.
  ############################################################
  if(is.null(n_modules)) {
    # Generate extra modules sizes. Should not need more than p.
    module_list <- vector("list", p)
  } else {
    module_list <- vector("list", n_modules)
  }
  
  all_nodes <- 1:p
  nodes_available <- min(2 * module_sizes[1], p)
  prob <- rep(1/p, p) # Initial probability of selecting a node for a module.
  deg <- rep(0, p)
  node_unselected <- rep(TRUE, p)
  need_more_modules <- TRUE
  i <- 0
  while(need_more_modules && i < p) {
    i <- i + 1
    m <- module_sizes[i] 
    index_nodes_available <- 1:nodes_available
    index_nodes_selected <- which(!node_unselected[index_nodes_available])
    nodes <- rep(0, m)
    if(i > 1) {
      link_node <- sample(index_nodes_selected, 1, 
                          prob = prob[index_nodes_selected])
      nodes[-m] <- sample(index_nodes_available[-link_node], m - 1, 
                          prob = prob[index_nodes_available][-link_node])
      nodes[m] <- link_node
    } else {
      nodes <- sample(index_nodes_available, m, prob = prob[index_nodes_available])
    }
    nodes <- sort(nodes)
    module_list[[i]] <- random_module(nodes, 
                                      deg_ex = deg[nodes], 
                                      ...)
    node_unselected[nodes] <- FALSE
    
    # Update degree distribution
    deg[nodes] <- deg[nodes] + apply(get_adjacency_matrix(module_list[[i]]), 2, sum)
    prob[!node_unselected] <- 1 / p * ecdf(deg[!node_unselected])(deg[!node_unselected])^selection_weight
    
    if(is.null(n_modules)) {
      need_more_modules <- any(node_unselected)
    } else {
      need_more_modules <- (i < n_modules)
    }
    nodes_available <- min(nodes_available + module_sizes[i], p)
  }
  
  # i = n_modules
  module_list <- module_list[1:i] # Remove any extra nodes.
  return(module_list)
}

###########################################################################
#
# Getter functions for network objects.
#
###########################################################################

#' Get the adjacency matrix of a network
#' 
#' The adjacency matrix indicates direct connections in the network. The ij'th
#' entry in the adjacency matrix is 1 if gene i and j are connected in any 
#' module and is 0 otherwise.
#' @param network A 'network' object; can be either weighted or unweighted.
#' @return An adjacency matrix with entry ij = 1 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero. 
#' @note Note that the connections in an adjacency matrix and association matrix
#' will likely differ. The adjacency matrix only considers direct connections in 
#' the network, whereas the association matrix takes into account the fact that 
#' there are overlapping modules which can create conditional dependencies
#' between two genes in seperate modules (i.e. genes that don't have a direct
#' connection in the graph).
#' @export
get_adjacency_matrix.network <- function(network, ...) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
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
  
  
  colnames(adj_matrix) <- network$node_names
  
  return(adj_matrix)
}


#' Get the association matrix of a network
#' 
#' The association matrix for a weighted network is recovered from its
#' covariance structure (see ?get_sigma.network). The off-diagonal entries are 
#' interpreted as partial correlations, and the diagonal entries are set to zero.
#' @param network A weighted 'network' object.
#' @return An association matrix with entry ij != 0 if node i and j have a nonzero
#' conditional linear association (in the Gaussian graphical model), and 0 
#' otherwise. The diagonal entries are all zero.
#' @note Note that the connections in an adjacency matrix and association matrix
#' will likely differ. The adjacency matrix only considers direct connections in 
#' the network, whereas the association matrix takes into account the fact that 
#' there are overlapping modules which can create conditional dependencies
#' between two genes in seperate modules (i.e. genes that don't have a direct
#' connection in the graph).
#' @export
get_association_matrix.network <- function(network, ...) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  if(!all(sapply(network$modules, is_weighted))) 
    stop("Network modules must be weighted.")
  
  sigma <- get_sigma(network)
  assoc_matrix <- -solve(sigma)
  diag(assoc_matrix) <- -diag(assoc_matrix)
  assoc_matrix <- cov2cor(assoc_matrix)
  
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
    stop("Network modules must be weighted.")
  
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
#' @return A list containing summary information for the node; this includes 
#' a vector of indicies to other nodes in the network it is connected to, and 
#' a vector of incidices to modules that contain the node.
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
  connections <- NULL
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
        module_index <- c(module_index, i)
        node_index_in_module <- which(network$modules[[i]]$nodes == node_index)
        adj_matrix <- get_adjacency_matrix(network$modules[[i]])
        connected_nodes <- network$modules[[i]]$nodes[
          which(adj_matrix[, node_index_in_module] != 0)]
        connections <- c(connections, connected_nodes)
      }
    }
  }
  
  # Delete any duplicated connections (this happens if two genes are connected
  # in multiple modules). 
  connections <- unique(connections)
  connections <- sort(connections)
  
  return(list(node = node,
              node_index = node_index,
              connections = connections,
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
  # Only perform checks if module_list is not NULL and non-empty.
  if(!is.null(module_list) && length(module_list) > 0) {
    if(class(module_list) == "network_module") {
      # A single module is provided; put into a list without warning.
      module_list <- list(module_list)
    }
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

#' Rewire connections to a node.
#' 
#' @param network A 'network' object to modify. 
#' @param node The node to rewire.
#' @param prob_rewire A value between 0 and 1. Each connection to 'node' 
#' will be rewired with probability equal to 'prob_rewire'. Note, the degree of 
#' 'node' is unchanged after this operation.
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling nodes to rewire.
#' @param exponent The exponent used for weighted sampling. When exponent = 0,
#' nodes are sampled uniformly. When exponent > 0, the sampling probability
#' is based on node weights.
#' @return The modified network.
#' @export
rewire_connections_to_node.network <- function(network,
                                               node,
                                               prob_rewire,
                                               weights = NULL,
                                               exponent = 0) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  if(length(network$modules) > 0) {
    for(i in 1:length(network$modules)) {
      if(node %in% network$modules[[i]]$nodes) {
        network$modules[[i]] <- 
          rewire_connections_to_node(network$modules[[i]], node,
                                     prob_rewire, weights, exponent)
      }
    }
  } else {
    warning("Argument 'network' contains no modules. Returning network unmodified.")
  }
  
  return(network)
}



#' Remove connections to a node.
#' 
#' @param network A 'network' object to modify. 
#' @param node The node to unwire.
#' @param prob_remove A value between 0 and 1. Each connection to 'node_index' 
#' will be removed with probability equal to 'prob_remove'.
#' @param weights (Optional) A vector of weights for each node. These are used
#' in addition to the degree of each node when sampling neighbors to unwire from.
#' @param exponent The exponent used for weighted sampling. When exponent = 0,
#' neighboring nodes are sampled uniformly. When exponent > 0, the sampling 
#' probability is based on node weights.
#' @return The modified network.
#' @export
remove_connections_to_node.network <- function(network,
                                               node,
                                               prob_remove,
                                               weights = NULL,
                                               exponent = 0) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  if(length(network$modules) > 0) {
    for(i in 1:length(network$modules)) {
      if(node %in%  network$modules[[i]]$nodes) {
        network$modules[[i]] <- 
          remove_connections_to_node(network$modules[[i]], node,
                                     prob_remove, weights, exponent)
      }
    }
  } else {
    warning("Argument 'network' contains no modules. Returning network unmodified.")
  }
  
  return(network)
}


#' Perturbs the connections in a network
#' 
#' The network is perturbed by removing connections from hubs and/or rewiring
#' other nodes in the network. By default, one hub is turned off (i.e. its 
#' connections are removed each with probability 'rewire_hub_prob' = 0.5), and 
#' no other nodes are changed. Hub nodes are defined as those having degree
#' above three standard deviations from the average degree, and nodes are
#' sampled from these to be turned off; if there are no hub nodes, then
#' those with the largest degree are turned off.
#' @param network The network to modify.
#' @param n_hubs The number of hub nodes to turn off.
#' @param n_nodes The number of non-hub nodes to rewire. When rewiring, the
#' degree of the node is unchanged.
#' @param rewire_hub_prob The probability that a connection is removed from
#' a hub that is selected to be turned off. If 'rewire_hub_prob' = 1, then
#' all of the connections to the hub are removed.
#' @param rewire_other_prob The probability that a connection is rewired from
#' a non-hub that is selected for rewiring. If 'rewire_other_prob' = 1, then 
#' all of the connections to the hub are rewired; however, this does not mean 
#' that all connections will be changed, as some connections may be removed
#' but later rewired back.
#' @param ... Additional arguments passed to rewire_connections_to_node() and 
#' remove_connections_to_node()
#' @return The modified network.
#' @export
perturb_network <- function(network, 
                            n_hubs = 1,
                            n_nodes = 0,
                            rewire_hub_prob = 0.5,
                            rewire_other_prob = 0.1,
                            ...) {
  if(!(class(network) == "network")) 
    stop(paste0("'", deparse(substitute(network)), 
                "' is not a 'network' object."))
  
  if(n_hubs > network$p || n_hubs < 0) {
    stop("Argument 'n_hubs' must be between 0 and the network size, p.")
  }
  if(n_nodes > network$p || n_nodes < 0) {
    stop("Argument 'n_nodes' must be between 0 and the network size, p.")
  }
  
  # First rewire others, then turn off the hub.
  # Hub genes will be identified by having a degree above 3 SDs from the mean.
  degrees <- get_degree_distribution(network)
  index_hubs <- which(degrees > floor(mean(degrees) + 3 * sd(degrees)))
  index_others <- setdiff(1:network$p, index_hubs)
  
  if(n_hubs > 0) {
    if(length(index_hubs) == 0) {
      # If no hubs, use nodes with highest degree.
      selected_hubs <- order(degrees, decreasing = TRUE)[1:n_hubs]
    } else if(length(index_hubs) < n_hubs) {
      # If not enough hubs, use hubs and nodes with highest degree.
      selected_hubs <- c(index_hubs, 
                         order(degrees, decreasing = TRUE)[1:(n_hubs - length(index_hubs))])
    } else {
      # Otherwise, sample from hubs.
      if(length(index_hubs) == 1) { # && n_hubs == 1, but don't need to check.
        # Do not use sample() if only 1 hub.
        selected_hubs <- index_hubs
      } else {
        selected_hubs <- sample(index_hubs, n_hubs)
      }
    }
    # Turn off selected hubs.
    for(index in selected_hubs) {
      network <- remove_connections_to_node(network, index, rewire_hub_prob, ...)
    }
  }
  
  if(n_nodes > 0) {
    if(length(index_others) < n_nodes) {
      # If not enough nodes, give warning and reduce n_nodes.
      warning(paste0(n_nodes, "nodes requested to rewire but only", length(index_others),
                     "were available and rewired."))
      n_nodes <- length(index_others)
    } 
    # Select nodes to rewire.
    if(length(index_others) == 1) {
      # Do not use sample() if only 1 hub.
      selected_nodes <- index_others
    } else {
      selected_nodes <- sample(index_others, n_nodes)
    }
    # Rewire selected nodes.
    for(index in selected_nodes) {
      network <- rewire_connections_to_node(network, index, rewire_other_prob, ...)
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
  
  n_modules <- length(network$modules)
  vals <- get_network_characteristics(network)
  
  message <- paste0(ifelse(is_weighted(network), "A weighted", "An unweighted"), 
                    " network containing ", vals$p, " nodes, ", 
                    vals$n_edges, " edges, and ", n_modules, 
                    " module", ifelse(n_modules > 1, "s", ""), ".\n")
  cat(message)
  print(round(unlist(vals[-c(1, 2, 6)]), 3))
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


#' Characteristics of the network topology
#' 
#' The average degree, clustering coefficient, and average path length are calculated.
#' @param network A 'network', 'network_module', or 'matrix' object.
#' @return A list containing characteristics of the network.
#' @export
get_network_characteristics <- function(network) {
  adj <- get_adjacency_matrix(network)
  graph <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  
  # Calculate the degree distribution P(K).
  deg <- igraph::degree(graph)
  tab <- table(deg)
  
  # Calculate the average clustering coefficient C(K) for nodes with each degree K.
  ind <- which(deg > 1)
  if(length(ind) == 0) {
    # Store properties of K in a data.frame.
    df <- data.frame(K = as.numeric(names(tab)), 
                     P_K = as.numeric(tab / sum(tab)),
                     C_K = NA)
    
  } else {
    clustco <- igraph::transitivity(graph, type = "local")[ind]
    deg_sub <- deg[ind]
    clus <- rep(NA, max(deg_sub))
    for(i in 2:(max(deg_sub))) {
      coefs <- clustco[which(deg_sub == i)]
      if(length(coefs) > 0) {
        clus[i] <- mean(coefs)
      }
    }
    
    # Store properties of K in a data.frame.
    if("1" %in% names(tab)) {
      C_K <- c(NA, clus[which(!is.na(clus))])
    } else {
      C_K <- clus[which(!is.na(clus))]
    }
    df <- data.frame(K = as.numeric(names(tab)), 
                     P_K = as.numeric(tab / sum(tab)),
                     C_K = C_K)
  }
  
  return(list(p = ncol(adj),
              n_edges = sum(adj[lower.tri(adj)]),
              `avg degree` = mean(deg), 
              `clustering coef` = igraph::transitivity(graph), 
              `avg path length` = igraph::mean_distance(graph), 
              df = df))
}
