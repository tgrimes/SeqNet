#' The Zero-Inflated Negative Binomial Distribution
#' 
#' @param x A vector of quantities.
#' @param size The dispersion paramater used in dnbinom.
#' @param mu The distribution mean.
#' @param rho The zero-inflation parameter.
#' @param log Logical; if TRUE, then log(d) is returned.
#' @param return The value(s) of the density function evaluated at x.
#' @export 
dzinb <- function(x, size, mu, rho = 0, log = FALSE) {
  d <- rho * (x == 0) + (1 - rho) * dnbinom(x, size, mu = mu)
  if(log) {
    d <- log(d)
  }
  return(d)
}

#' The Zero-Inflated Negative Binomial Distribution
#'
#' @param p A vector of probabilities
#' @param size The dispersion paramater used in dnbinom.
#' @param mu The distribution mean.
#' @param rho The zero-inflation parameter.
#' @param lower.tail Logical; if TRUE, then probabilities are P(X <= x).
#' Otherwise, P(X > x).
#' @param log.p Logical; if TRUE, then exp(p) is used.
#' @export 
qzinb <- function (p, size, mu, rho, lower.tail = TRUE, log.p = FALSE) {
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- exp(p)
  }
  x <- qnbinom(pmax(0, (p - rho)/(1 - rho)), size, mu = mu)
  
  return(x)
}

#' The Zero-Inflated Negative Binomial Distribution
#'
#' @param q A vector of quantities.
#' @param size The dispersion paramater used in dnbinom.
#' @param mu The distribution mean.
#' @param rho The zero-inflation parameter.
#' @param lower.tail Logical; if TRUE, then probabilities are P(X <= x).
#' Otherwise, P(X > x).
#' @param log.p Logical; if TRUE, then log(p) is returned.
#' @export 
pzinb <- function (q, size, mu, rho, lower.tail = TRUE, log.p = FALSE) {
  p <- rho * (q >= 0) + (1 - rho) * pnbinom(q, size, mu = mu)
  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  
  return(p)
}

#' The Zero-Inflated Negative Binomial Distribution
#'
#' @param n The number of random values to return.
#' @param size The dispersion paramater used in dnbinom.
#' @param mu The distribution mean.
#' @param rho The zero-inflation parameter.
#' @export 
rzinb <- function (n, size, mu, rho) {
  x <- ifelse(rbinom(n, 1, rho), 0, rnbinom(n, size, mu = mu))
  
  return(x)
}

#' Generate RNA-seq counts
#' 
#' The count data are generated based on the gene-gene associations of an
#' udnerlying network. An association structure is imposed by first generating 
#' data from a multivariate Gaussian distribution, and counts are then obtained
#' through the inverse tranformation method. To generate realistic counts, either 
#' a reference dataset or parameters for the ZINB model (size, mu, rho) can be provided.
#' @param n The number of samples to generate.
#' @param network A 'network' object.
#' @param reference Either a vector or data.frame of counts from a reference
#' gene expression profile. If a data.frame is provided, each column should
#' correspond to a gene. If both 'reference' and 'params' are NULL, then parameters
#' are estimated from the kidney dataset.
#' @param params A matrix of ZINB parameter values; each column should contain 
#' the size, mu, and rho parameters for a gene. 
#' @param verbose Boolean indicator for message output.
#' @return A list containing the generated counts and the ZINB parameters used
#' to create them.
#' @export 
gen_counts <- function(n, 
                       ...,
                       reference = NULL,
                       params = NULL,
                       verbose = TRUE) {
  if(n <= 0) {
    stop("n must be positive.")
  }
  
  if(is.null(reference) && is.null(params)) {
    warning("Using kidney data as reference dataset.")
    reference <- get_kidney_reference_data()
    reference <- sample_reference_data(reference, p)
  }
  
  if(!is.null(params) && network_list[[1]]$p != ncol(params)) {
    stop("'params' must have the same number of columns as genes in the network(s).")
  }

  # Estimate model paramters from reference dataset
  if(is.null(params)) {
    index <- 1:p # Default subset for columns in reference.
    if(p > ncol(reference)) {
      if(verbose) {
        cat("reference contains fewer columns than nodes in network.",
            "Sampling genes with replacement.\n")
      }
      index <- sample(1:ncol(reference), p, replace = TRUE)
    } else if (p < ncol(reference)) {
      if(verbose) {
        cat("reference contains more columns than nodes in network.",
            "Sampling a subset of genes.\n")
      }
      index <- sample(1:ncol(reference), p)
    }
    # Subset the columns of the reference.
    reference <- reference[, index]
    model <- est_params_from_reference(
      reference,
      verbose = verbose)
    params <- model$params
  } # End if(is.null(params))
  
  colnames(params) <- colnames(reference)
  
  gen <- gen_gaussian(n, network)
  x <- pnorm(gen$x) # Obtain n by p matrix of quantiles.
  
  # Convert quantiles to counts for each gene.
  for(i in 1:p) {
    x[, i] <- qzinb(x[, i], 
                    size = params[1, i], 
                    mu = params[2, i],
                    rho = params[3, i])
  }
  
  return(list(x = x,
              params = params))
}




#' Estimate ZINB parameters from reference data
#' 
#' The observations in the reference dataset should be as homogeneous as possible.
#' For example, we should not expect differential expression or differential 
#' connectivity of genes within the sample. If the data are heterogeneous, the
#' estimation of the parameters may be unreliable.
#' @param reference Either a vector or data.frame of counts from a reference
#' gene expression profile. If a data.frame is provided, each column should
#' correspond to a gene. If NULL, the kidney dataset is used and parameter estimates
#' for all of its genes are calculated.
#' @param verbose Boolean indicator for message output.
#' @return Returns a list containing a matrix of parameter estimates 'size', 
#' 'mu', and 'rho' for each gene in the reference, and the reference dataset
#' used.
#' @export
est_params_from_reference <- function(reference = NULL,
                                      verbose = TRUE) {
  get_nb_params <- function(x) {
    fit <- fitdistrplus::fitdist(x, "nbinom",
                                 lower = 0)$estimate
    params <- list(size = unname(fit)[1], 
                   mu = unname(fit)[2],
                   rho = 0,
                   BIC = NA)
    L <- dnbinom(x, size = params$size, mu = params$mu)
    bic_nb <- -2 * sum(log(L)) + log(length(x)) * 2 # Two parameters in model.
    params$BIC <- bic_nb
    
    return(params)
  }
  
  get_zinb_params <- function(x) {
    # Use estimated size and mu from nb model to start zinb.
    fitzi <- NULL
    perc_zero <- mean(x == 0)
    if(perc_zero > 0) {
      # If x contains any zero values, try fitting a ZI model.
      fitzi <- tryCatch({
          m <- mean(x[x != 0]) # Used to calculate initial value.
          fit <- fitdistrplus::fitdist(x, dzinb, 
                                       start = list(size = m^2 / (var(x) - m), 
                                                    mu = m, 
                                                    rho = perc_zero / 2),
                                       lower = 0)$estimate
          if(fit[3] <= 0) {
            # If estimate for rho (ZI parameter) is not positive, return NULL.
            fit <- NULL
          } 
          fit
        }, error = function(e) {
          warning(paste("Estimation of ZINB parameters failed:", e))
          return(NULL)
        })
    }
    if(is.null(fitzi)) {
      # If estimation of ZINB failed or was not performed, return with BIC at Inf.
      params <- list(size = NA, mu = NA, rho = NA, BIC = Inf)
    } else {
      params <- list(size = unname(fitzi[1]), 
                     mu = unname(fitzi[2]),
                     rho = unname(fitzi[3]),
                     BIC = NA)
      L <- dzinb(x, params$size, params$mu, params$rho)
      params$BIC <- -2 * sum(log(L)) + log(length(x)) * 3 # Three parameters in model.
    }
    
    return(params)
  }
  
  get_params <- function(x) {
    # First fit NB model.
    fit <- get_nb_params(x)
    if(sum(x == 0) > 0) {
      # If there are zeroes in data, compare NB to ZINB using BIC criterion.
      fit_zinb <- get_zinb_params(x)
      if((fit_zinb$BIC <= fit$BIC) && (fit_zinb$rho > 0)) {
        # If ZINB is better fit and estimated rho is positive, use ZINB model.
        fit <- fit_zinb
      }
    }
    
    # # Get shape and scale parameters for gamma-poisson.
    # est_prob <- fit$size / (fit$size + fit$mu)
    # est_scale <- (1 - est_prob) / est_prob
    
    params <- c(size = fit$size,
                mu = fit$mu,
                rho = fit$rho)
    
    return(params)
  }
  
  if(is.null(reference)) {
    if(verbose) {
      cat("No reference dataset provided. Using kidney data.\n")
    }
    reference <- get_kidney_reference_data()
  }
  
  if(is.vector(reference)) {
    reference <- matrix(reference, ncol = 1)
  } else if(!is.data.frame(reference) && !is.matrix(reference)) {
    stop("reference should be a vector or data.frame.")
  }
  
  p <- ncol(reference)
  params <- matrix(0, nrow = 3, ncol = p)
  colnames(params) <- colnames(reference)
  rownames(params) <- c("size", "mu", "rho")
  
  for(i in 1:p) {
    x <- reference[, i]
    params[, i] <- get_params(x)
  }
  
  return(list(params = params,
              reference = reference))
}


#' Get reference data from kidney dataset
#' 
#' @return A data.frame containing a reference dataset with each column 
#' containing an observed expression profile for a gene.
#' @export 
get_kidney_reference_data <- function() {
  data("kidney")
  # reference <- t(kidney$counts)
  # Only use samples from non-tumor tissue. This avoids differential expression.
  kidney <- kidney$tumor
  
  # Remove any genes that have zero variation or have extremely low expression.
  index <- apply(kidney, 2, function(x) sd(x) > 0)
  kidney <- kidney[, index]
  index <- apply(kidney, 2, function(x) mean(x) > 5)
  kidney <- kidney[, index]
  
  return(kidney)
}


#' Sample genes from reference dataset
#' 
#' @param reference_data The reference data.frame to use.
#' @param p The number of genes (columns) to sample
#' @param percent_ZI The percentage of genes to be zero inflated. If
#' NULL, the genes are sampled at random; in this case, the empirical 
#' distribution of gene expression profiles will determine the probablility
#' that a sampled gene is zero inflated.
#' @param min_mean Genes with mean expression below this value are removed.
#' @param max_mean_percentile A value between 0 and 1. All genes with mean 
#' expression above this percentile are removed.
#' @param threshold_ZI The minimum proportion of zero counts for a gene to be
#' considered as zero inflated.
#' @export 
sample_reference_data <- function(reference_data,
                                  p,
                                  percent_ZI = NULL,
                                  min_mean = 30,
                                  max_mean_percentile = 0.99,
                                  threshold_ZI = 0.2) {
  if(p <= 0) 
    stop("p must be greater than 0.")
  if(min_mean < 0) 
    stop("min_mean must be greater than 0.")
  if(max_mean_percentile < 0 || max_mean_percentile > 1)
    stop("max_mean_percentile must be between 0 and 1.")
  
  means <- apply(reference_data, 2, mean)
  max_mean <- quantile(means, max_mean_percentile)
  reference_data <- reference_data[, (means >= min_mean) & (means <= max_mean)]

  index_zero <- NULL
  index_ZI_genes <- NULL
  if(!is.null(percent_ZI)) {
    if(percent_ZI > 1 || percent_ZI < 0) {
      stop("percent_ZI must be between 0 and 1.")
    }
    # index_ZI_genes <- which(apply(reference_data, 2, 
    #                               function(x) mean(x == 0) > threshold_ZI))
    index_ZI_genes <- which(apply(reference_data, 2, function(x) {
      perc_zero <- mean(x == 0)
      if(perc_zero < 0.05) {
        return(FALSE)
      } else if(mean(x) < 30) {
        return(FALSE) 
      } else if(var(x) == 0){
        return(FALSE)
      }
      tryCatch(
        {
          m <- mean(x[x != 0]) # Used to calculate initial value.
          fit_zinb <- fitdistrplus::fitdist(x, dzinb, 
                                            start = list(size = m^2 / (var(x) - m), 
                                                         mu = m, 
                                                         rho = perc_zero / 2),
                                            lower = c(1, 0, 0))
          fit_nb <- fitdistrplus::fitdist(x, dzinb, 
                                          start = list(size = m^2 / (var(x) - m), 
                                                       mu = m),
                                          lower = c(0, 0))
          L <- dzinb(x, fit_zinb$estimate[1], fit_zinb$estimate[2], fit_zinb$estimate[3])
          bic_zinb <- -2 * sum(log(L)) + log(length(x)) * 3 # Three parameters in model.
          L <- dzinb(x, fit_nb$estimate[1], fit_nb$estimate[2])
          bic_nb <- -2 * sum(log(L)) + log(length(x)) * 2 # Two parameters in model.
          return(bic_zinb < bic_nb)
        },
        error = function(e) {
          return(FALSE)
        })
    }))
    index_ZI_genes <- unname(index_ZI_genes)
    p_0 <- ceiling(p * percent_ZI)
    if(p_0 > 0) {
      if(length(index_ZI_genes) == 0) {
        stop("None of the genes in the reference dataset are zero inflated.")
      } else if(length(index_ZI_genes) < p_0) {
        warning(paste("Not enough zero-inflated genes in reference dataset.",
                      "Sampling with replacement will be used."))
        index_zero <- sample(index_ZI_genes, p_0, replace = TRUE)
      } else {
        index_zero <- sample(index_ZI_genes, p_0)
      }
    }
    p <- p - p_0
  }
  
  index_nonzero <- NULL
  if(p > 0) {
    # If percent_zero_inflated is NULL, then index_genes will contain all genes.
    # Otherwise, the zero inflated genes will be removed.
    index_genes <- setdiff(1:ncol(reference_data), index_ZI_genes)
    # Note, p will be zero if ceiling(p * percent_zero_inflated) is 1.
    if(p > 0)
      if(length(index_genes) == 0) {
        stop("None of the genes in the reference dataset are non-zero-inflated.")
      } else if(length(index_genes) < p) {
        warning(paste("Not enough genes in reference dataset.",
                      "Sampling with replacement will be used."))
        index_nonzero <- sample(index_genes, p, replace = TRUE)
      } else {
        index_nonzero <- sample(index_genes, p)
      }
  }
  
  index <- c(index_zero, index_nonzero)
  
  return(reference_data[, index])
}
