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
get_adjacency_matrix.default <- function(...) {
  cat("get_adjacency_matrix() is defined for 'network' and 'network_module' objects.\n")
}
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
#' The association matrix is constructed from all modules in a network.
#' @param ... Either a 'network', 'network_module', or 'matrix' object.
#' @return An association matrix with entry ij != 0 if node i and j are 
#' connected, and 0 otherwise. The diagonal entries are all zero.
#' @export
get_association_matrix <- function(...) {
  UseMethod("get_association_matrix")
}
get_association_matrix.default <- function(...) {
  cat("get_association_matrix() is defined for 'network' and 'network_module' objects.\n")
}
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
  
  # Set any values that are near zero to exactly zero.
  m[abs(m) <= 10^-13] <- 0
  
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
get_sigma.default <- function(...) {
  cat("get_sigma() is defined for 'network' and 'network_module' objects.\n")
}
get_sigma.matrix <- function(m, ...) {
  if(any(diag(m) == 0)) {
    stop("Matrix 'm' should have nonzero entries along its diagonal.")
  }
  warning("Interpreting matrix 'm' as a partial correlation matrix.")
  assoc_matrix <- m
  assoc_matrix <- -assoc_matrix
  diag(assoc_matrix) <- -diag(assoc_matrix)
  sigma <- solve(assoc_matrix)
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
is_weighted.default <- function(...) {
  cat("is_weighted() is defined for 'network' and 'network_module' objects.\n")
}
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
  
  if(!all(m %in% c(0, 1))) {
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
remove_weights.default <- function(...) {
  cat("remove_weights() is defined for 'network', 'network_module', and 'matrix' objects.\n")
}
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
get_node_names.default <- function(...) {
  cat("get_node_names() is defined for 'network' and 'network_module' objects.\n")
}
get_node_names.matrix <- function(m, ...) {
  return(colnames(m))
}