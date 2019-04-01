# setGeneric(
#   name = "get_adjacency_matrix", 
#   def = function(object) standardGeneric("get_adjacency_matrix")
# )
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

# setGeneric(
#   name = "get_association_matrix", 
#   def = function(object) standardGeneric("get_association_matrix")
# )
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


remove_weights <- function(...) {
  UseMethod("remove_weights")
}
remove_weights.default <- function(...) {
  cat("remove_weights() is defined for 'network' and 'network_module' objects.\n")
}


get_node_names <- function(...) {
  UseMethod("get_node_names")
}
get_node_names.default <- function(...) {
  cat("get_node_names() is defined for 'network' and 'network_module' objects.\n")
}
get_node_names.matrix <- function(m, ...) {
  return(colnames(m))
}

cov2prec <- function(m) {
  prec <- -solve(m)
  diag(prec) <- -diag(prec)
}