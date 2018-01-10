
# library(doParallel)
# library(foreach)
# library(mgcv) # gam() for fitting nb
# library(Matrix) # crossprod()

#' Combined partial least squares with negative-binomial
#' 
#' Reconstructs a full genomic association network from a set of 
#' next-generation sequencing count data. No preprocessing is performed on x. 
#' @param x The n by p matrix of raw count data. Rows should correspond to
#'   samples, and columns should correspond to genes.
#' @param v The number of PLS components to compute.
#' @param threshold (optional) the fdr signifiance threshold for inferring the 
#'   network through empirical bayes fdr. If NULL, all scores are returned. 
#'   Otherwise, scores with fdr above this threshold are set to zero.
#' @param parallel Should cpls be run in parallel? (Unix system required.)
#' @return A list containing the p by p matrix of association scores, and the 
#' threshold used to determine significant associations.
#' @export
run_cpls <- function(x, threshold = NULL, v = 3, parallel = FALSE) {
  if(parallel) {
    if(.Platform$OS.type == "unix") {
      return(run_cpls_parallel(x = x, v = v, threshold = threshold))
    } else {
      cat("OS must be Unix to run in parallel. Running `cpls` sequentially.\n")
    }
  }
  
  n <- nrow(x) #Number of observations.
  p <- ncol(x) #Number of genes.
  if(v >= n) {
    stop("v must be less than n (the number of rows in X)")
  }
  if(n > p) {
    warning("expected p > n; matrix should have observations in rows and 
            genes in columns.")
  }
  if(is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if(any(!is.numeric(x))) {
    stop("x contains non-numeric value(s).")
  }
  
  #Total read counts in each sample.
  num_reads_by_sample <- apply(x, 1, sum) 
  
  #The association scores.
  s <- matrix(0, nrow = p, ncol = p)
  
  #Standardization of explanatory variables.
  standardize <- function(x) {
    apply(x, 2, function(x) { (x - mean(x))/sd(x) })
  }
  #If any gene has zero variance, add one count to its first sample.
  zero_variance_genes <- which(apply(x, 2, function(x) var(x) == 0))
  if(length(zero_variance_genes) > 0) {
    cat("Adding one count to first sample for genes:", zero_variance_genes, ".")
    x[1, zero_variance_genes] <- x[1, zero_variance_genes] + 1
  }
  x_std <- standardize(x)
  
  
  #S1-S5.
  for(i in 1:p) {
    if((i %% 1000) == 0) cat("gene", i, "of", p, "\n")
    
    #The PLS coefficients and components.
    pls_coef <- array(0, dim = c(p, v)) 
    pls_comp <- array(0, dim = c(n, v))
    
    #Initialize the deflated x matrix used in each iteration.
    #These explanatory variables ARE standardized.
    x_deflated <- x_std[, -i]
    
    #S2: compute the PLS coefficients and components.
    for(j in 1:v) {
      #The response variable IS standardized.
      pls_coef[-i, j] <- Matrix::crossprod(x_deflated, x[, i]) # t(x_deflated) %*% x[, i]
      pls_coef[, j]  <- pls_coef[, j]  / sqrt(sum(pls_coef[, j]^2))
      pls_comp[, j] <- x_deflated %*% pls_coef[-i, j]
      
      temp <- tryCatch(solve(Matrix::crossprod(pls_comp[, j])), 
                       error = function(e) { warning(e); return(NA) })
      #Check that solve() was successful.
      if(!is.na(temp)) {
        x_deflated <- x_deflated - 
          pls_comp[, j] %*% temp %*% Matrix::crossprod(pls_comp[, j], x_deflated)
      } else {
        #Only stop if fewer than v components have been computed.
        if(j < v) {
          stop(paste("cPLS algorithm halted. Only v =", j, "of", 
                     v, "components computed."))
        }
      }
    }
    
    #S3: fit a negative binomial model: log(E(X)) = pls_comp'beta + log(N).
    #The negative binomial model coefficients.
    beta <- array(0, dim = c(v))
    
    #The response variable is NOT standardized.
    temp_data <- data.frame(y = x[, i], 
                            C = pls_comp, 
                            total = num_reads_by_sample)
    form <- as.formula(paste("y ~", paste(colnames(temp_data[2:(v + 1)]), 
                                          collapse = "+")))
    
    beta <- tryCatch(
      coef(mgcv::gam(form, offset = log(total), family = mgcv::nb(), data = temp_data))[-1],
      error = function(e) {
        print(paste(e, "- gene", i, "of", p, "- setting coefs to zero."))
        return(rep(0, v))
      })
    
    #S4: compute the association scores for each gene pair i and k,
    #    k from 1 to p.
    s[, i] <- pls_coef %*% beta
  }
  
  #S6: Set the diagonal of s to 1 and symmetrize s.
  diag(s) <- 0
  s <- (s + t(s)) * 0.5
  
  # If threshold is provided, set scores with fdr above this threshold to zero.
  if(!is.null(threshold)) {
    cat("Computing association matrix using fdr rate", threshold, "\n")
    s[fdr(normalize_cpls_scores(s))$likelihood > threshold] <- 0
  }
  
  return(list(scores = s, threshold = threshold))
}

#' Normalize the cPLS scores
#' 
#' Performs the transformation: sign(x) log(|x| + 1). This normalizes the
#' cPLS scores; should be performed prior to using `fdr()`.
#' @param scores The association scores from the cPLS algorithm.
#' @export
normalize_cpls_scores <- function(scores) {
  sign(scores) * log(abs(scores) + 1)
}

#' 
#' 
#' #' Combined partial least squares with negative-binomial
#' #' 
#' #' Parallel version of run_cpls(). Reconstructs a full genomic association 
#' #' network from a set of next-generation sequencing count data.
#' #' @param x The n by p matrix of raw count data. Rows should correspond to
#' #'   samples, and columns should correspond to genes.
#' #' @param v The number of PLS components to compute.
#' #' @param threshold (optional) the fdr signifiance threshold for inferring the 
#' #'   network through empirical bayes fdr. If NULL, all scores are returned. 
#' #'   Otherwise, scores with fdr above this threshold are set to zero.
#' #' @param parallel Should cpls be run in parallel? (Unix system required.)
#' #' @return a matrix of raw association scores. 
#' #' @export
#' run_cpls_parallel <- function(x, v = 3, threshold = NULL) {
#'   n <- nrow(x) #Number of observations.
#'   p <- ncol(x) #Number of genes.
#'   if(v >= n) {
#'     stop("v must be less than n (the number of rows in X)")
#'   }
#'   if(n > p) {
#'     warning("expected p > n; matrix should have observations in rows and 
#'             genes in columns.")
#'   }
#'   if(is.data.frame(x)) {
#'     x <- as.matrix(x)
#'   }
#'   if(any(!is.numeric(x))) {
#'     stop("x contains non-numeric value(s).")
#'   }
#'   
#'   #Total read counts in each sample.
#'   num_reads_by_sample <- apply(x, 1, sum) 
#'   
#'   #Standardization of explanatory variables.
#'   standardize <- function(x) {
#'     apply(x, 2, function(x) { (x - mean(x))/sd(x) })
#'   }
#'   #If any gene has zero variance, add one count to its first sample.
#'   zero_variance_genes <- which(apply(x, 2, function(x) var(x) == 0))
#'   if(length(zero_variance_genes) > 0) {
#'     cat("Adding one count to first sample for genes:", zero_variance_genes, ".")
#'     x[1, zero_variance_genes] <- x[1, zero_variance_genes] + 1
#'   }
#'   x_std <- standardize(x)
#'   
#'   #The association scores.
#'   s <- matrix(0, nrow = p, ncol = p)
#'   s_list <- vector("list", p)
#'   
#'   #S1-S5.
#'   n_cores <- parallel::detectCores()
#'   cat("Running cPLS on", n_cores, "cores.\n")
#'   doParallel::registerDoParallel(n_cores)
#'   s_list <- foreach::foreach(i = 1:p) %dopar% {
#'     
#'     #The PLS coefficients and components.
#'     pls_coef <- array(0, dim = c(p, v)) 
#'     pls_comp <- array(0, dim = c(n, v))
#'     
#'     #Initialize the deflated x matrix used in each iteration.
#'     #These explanatory variables ARE standardized.
#'     x_deflated <- x_std[, -i]
#'     
#'     #S2: compute the PLS coefficients and components.
#'     for(j in 1:v) {
#'       #The response variable IS standardized.
#'       pls_coef[-i, j] <- Matrix::crossprod(x_deflated, x[, i]) # t(x_deflated) %*% x[, i]
#'       pls_coef[, j]  <- pls_coef[, j]  / sqrt(sum(pls_coef[, j]^2))
#'       pls_comp[, j] <- x_deflated %*% pls_coef[-i, j]
#'       
#'       temp <- tryCatch(solve(Matrix::crossprod(pls_comp[, j])), 
#'                        error = function(e) { warning(e); return(NA) })
#'       #Check that solve() was successful.
#'       if(!is.na(temp)) {
#'         x_deflated <- x_deflated - 
#'           pls_comp[, j] %*% temp %*% Matrix::crossprod(pls_comp[, j], x_deflated)
#'       } else {
#'         #Only stop if fewer than v components have been computed.
#'         if(j < v) {
#'           stop(paste("cPLS algorithm halted. Only v =", j, "of", 
#'                      v, "components computed."))
#'         }
#'       }
#'     }
#'     
#'     #S3: fit a negative binomial model: log(E(X)) = pls_comp'beta + log(N).
#'     #The negative binomial model coefficients.
#'     beta <- array(0, dim = c(v))
#'     
#'     #The response variable is NOT standardized.
#'     temp_data <- data.frame(y = x[, i], 
#'                             C = pls_comp, 
#'                             total = num_reads_by_sample)
#'     form <- as.formula(paste("y ~", paste(colnames(temp_data[2:(v + 1)]), 
#'                                           collapse = "+")))
#'     
#'     beta <- tryCatch(
#'       coef(mgcv::gam(form, offset = log(total), family = mgcv::nb(), data = temp_data))[-1],
#'       error = function(e) {
#'         print(paste(e, "- gene", i, "of", p, "- setting coefs to zero."))
#'         return(rep(0, v))
#'       })
#'     
#'     #S4: compute the association scores for each gene pair i and k,
#'     #    k from 1 to p.
#'     return(pls_coef %*% beta)
#'   }
#'   doParallel::stopImplicitCluster()
#'   s <- sapply(s_list, cbind)
#'   
#'   #S6: Set the diagonal of s to 1 and symmetrize s.
#'   diag(s) <- 1
#'   s <- (s + t(s)) * 0.5
#'   
#'   # If threshold is provided, set scores with fdr above this threshold to zero.
#'   if(!is.null(threshold)) {
#'     cat("Computing association matrix using threshold rate", threshold, "\n")
#'     s[fdr(normalize_cpls_scores(s))$likelihood > threshold] <- 0
#'   }
#'   
#'   results <- list(scores = s, threshold = threshold)
#'   
#'   return(results)
#' }