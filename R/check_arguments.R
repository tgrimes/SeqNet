#' Wrapper for ArgumentCheck::finishArgCheck
#' 
#' @param check An object of class ArgCheck.
#' @details see ?ArgumentCheck::finishArgCheck for details. 
report_checks <- function(check) {
  ArgumentCheck::finishArgCheck(check)
}

#' Wrapper for ArgumentCheck::newArgCheck()
#' 
#' @details see ?ArgumentCheck::newArgCheck() for details. 
new_checklist <- function() {
  ArgumentCheck::newArgCheck()
}

#' Wrapper for ArgumentCheck::addError()
#' 
#' @details see ?ArgumentCheck::addError() for details. 
check_adjacency_matrix <- function(adjacency_matrix, checklist) {
  arg <- deparse(substitute(adjacency_matrix))
  
  if(!isSymmetric(unname(adjacency_matrix))) 
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must be symmetric."),
      argcheck = checklist
    )
  if(!all(adjacency_matrix %in% c(0, 1))) 
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' contains elements other than 0 or 1."),
      argcheck = checklist
    )
  if(any(is.na(adjacency_matrix))) 
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' contains NA values."),
      argcheck = checklist
    )
  
  return(checklist)
}

#' Wrapper for ArgumentCheck::addError()
#' 
#' @details see ?ArgumentCheck::addError() for details. 
check_association_matrix <- function(association_matrix, checklist) {
  arg <- deparse(substitute(association_matrix))
  
  if(!isSymmetric(unname(association_matrix))) 
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must be symmetric."),
      argcheck = checklist
    )
  if(any((abs(association_matrix) < 10^-14) && (association_matrix != 0))) 
    ArgumentCheck::addWarning(
      msg = paste0("Some values in argument '", arg, "' are near zero (< 10^-14). ",
                   "Consider setting these to zero."),
      argcheck = checklist
    )
  if(any(is.na(association_matrix))) 
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' contains NA values."),
      argcheck = checklist
    )
  
  return(checklist)
}

#' Wrapper for ArgumentCheck::addError()
#' 
#' @details see ?ArgumentCheck::addError() for details. 
check_precision_matrix <- function(precision_matrix, checklist) {
  arg <- deparse(substitute(precision_matrix))
  
  if(!isSymmetric(unname(precision_matrix))) 
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must be symmetric."),
      argcheck = checklist
    )
  if(any(abs(precision_matrix) < 10^-14 && precision_matrix != 0)) 
    ArgumentCheck::addWarning(
      msg = paste0("Some values in argument '", arg, "' are near zero (< 10^-14). ",
                  "Consider setting these to zero."),
      argcheck = checklist
    )
  if(!any(precision_matrix == 0)) 
    ArgumentCheck::addWarning(
      msg = paste0("Argument '", arg, "' contains all nonzero values.",
                  "Is covariance matrix is being used?"),
      argcheck = checklist
    )
  if(any(is.na(precision_matrix))) 
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' contains NA values."),
      argcheck = checklist
    )
  if(any(diag(precision_matrix) == 0)) 
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' contains zeroes along its diagonal."),
      argcheck = checklist
    )
  
  return(checklist)
}

#' Wrapper for ArgumentCheck::addError()
#' 
#' @details see ?ArgumentCheck::addError() for details. 
check_positive_integer <- function(p, checklist) {
  arg <- deparse(substitute(p))
  
  if(!is.numeric(p)) {
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must be a positive numeric value."),
      argcheck = checklist
    )
  } else if(p %% 1 != 0) {
      ArgumentCheck::addError(
        msg = paste0("Argument '", arg, "' = ", p, " must be an integer."),
        argcheck = checklist
      )
  }
  if(p <= 0) 
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must be positive."),
      argcheck = checklist
    )
  
  return(checklist)
}

#' Wrapper for ArgumentCheck::addError()
#' 
#' @details see ?ArgumentCheck::addError() for details. 
check_nonnegative_integer <- function(p, checklist) {
  arg <- deparse(substitute(p))
  
  if(!is.numeric(p)) {
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must be a positive numeric value."),
      argcheck = checklist
    )
  } else if(p %% 1 != 0) {
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' = ", p, " must be an integer."),
      argcheck = checklist
    )
  }
  if(p < 0) 
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must be non-negative"),
      argcheck = checklist
    )
  
  return(checklist)
}

#' Wrapper for ArgumentCheck::addError()
#' 
#' @details see ?ArgumentCheck::addError() for details. 
check_positive_integer_vector <- function(p, checklist) {
  arg <- deparse(substitute(p))
  
  if(!is.numeric(p)) {
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must contain positive numeric values."),
      argcheck = checklist
    )
  } else if(any(p %% 1 != 0)) {
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must contain integer values."),
      argcheck = checklist
    )
  }
  if(any(p <= 0))
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must contain all positive values."),
      argcheck = checklist
    )
  
  return(checklist)
}

#' Wrapper for ArgumentCheck::addError()
#' 
#' @details see ?ArgumentCheck::addError() for details. 
check_in_closed_interval <- function(x, checklist, 
                                      lower_bound, upper_bound) {
  arg <- deparse(substitute(x))
  if(x < lower_bound || x > upper_bound)
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must be in the interval [", 
                   lower_bound, ", ", upper_bound, "]."),
      argcheck = checklist
    )
}

#' Wrapper for ArgumentCheck::addError()
#' 
#' @details see ?ArgumentCheck::addError() for details. 
check_in_open_interval <- function(x, checklist, 
                                   lower_bound, upper_bound) {
  arg <- deparse(substitute(x))
  if(x <= lower_bound || x >= upper_bound)
    ArgumentCheck::addError(
      msg = paste0("Argument '", arg, "' must be in the interval (", 
                   lower_bound, ", ", upper_bound, ")."),
      argcheck = checklist
    )
}

