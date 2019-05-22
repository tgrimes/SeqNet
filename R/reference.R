
#' Get reference data from kidney dataset
#' 
#' @return A data.frame containing a reference dataset with each column 
#' containing an observed expression profile for a gene.
#' @export 
get_kidney_reference_data <- function() {
  reference <- readRDS("data/kidney.rds")
  
  # Remove any genes that have zero variation or have extremely low expression.
  has_variation <- apply(reference$expression, 2, function(x) length(unique(x)) != 1)
  reference$expression <- reference$expression[, has_variation]
  reference$gene_length <- reference$gene_length[has_variation]
  has_expression <- apply(reference$expression, 2, function(x) mean(x) > 5)
  reference$expression <- reference$expression[, has_expression]
  reference$gene_length <- reference$gene_length[has_expression]
  
  return(reference)
}

make_kidney_reference_data <- function() {
  data("kidney")
  # reference <- t(kidney$counts)
  # Only use samples from non-tumor tissue. This avoids differential expression.
  kidney <- kidney$normal
  
  # Remove any genes that have zero variation or have extremely low expression.
  has_variation <- apply(kidney, 2, function(x) length(unique(x)) != 1)
  kidney <- kidney[, has_variation]
  has_expression <- apply(kidney, 2, function(x) mean(x) > 5)
  kidney <- kidney[, has_expression]
  
  genes <- data.frame(entrezgene = colnames(kidney), stringsAsFactors = FALSE)
  genes <- entrez_to_lengths(genes, "human")
  saveRDS(genes, "data/genes.rds")
  
  index <- which(is.na(genes$length))
  # 1205 genes have NA lengths. Most of these appear to be "pseudogenes".
  
  if(length(index) > 0) {
    kidney <- kidney[, -index]
    genes <- genes[-index, ]
  }
 
  lengths <- genes$length
  names(lengths) <- genes$entrezgene
  
  reference <- list(expression = kidney, gene_length = length)
  saveRDS(reference, "data/kidney.rds")
  return(reference)
}


#' Sample genes from reference dataset
#' 
#' @param reference_data The reference data.frame to use.
#' @param p The number of genes (columns) to sample
#' @param percent_ZI The percentage of genes to be zero inflated. If
#' NULL, the genes are sampled at random; in this case, the empirical 
#' distribution of gene expression profiles will determine the probablility
#' that a sampled gene is zero inflated.
#' @param threshold_ZI The minimum proportion of zero counts for a gene to be
#' considered as zero inflated.
#' @export 
sample_reference_data <- function(reference_data,
                                  p,
                                  percent_ZI = NULL,
                                  threshold_ZI = 0.2) {
  if(p <= 0) 
    stop("p must be greater than 0.")
  
  index_zero <- NULL
  index_ZI_genes <- NULL
  if(!is.null(percent_ZI)) {
    if(percent_ZI > 1 || percent_ZI < 0) {
      stop("percent_ZI must be between 0 and 1.")
    }
    index_ZI_genes <- which(apply(reference_data, 2,
                                  function(x) mean(x == 0) > threshold_ZI))
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





# sample_reference_data <- function(reference_data,
#                                   p,
#                                   percent_ZI = NULL,
#                                   min_mean = 30,
#                                   max_mean_percentile = 0.99,
#                                   threshold_ZI = 0.2) {
#
# @param min_mean Genes with mean expression below this value are removed.
# @param max_mean_percentile A value between 0 and 1. All genes with mean 
# expression above this percentile are removed.
# 
# 
#   if(p <= 0) 
#     stop("p must be greater than 0.")
#   if(min_mean < 0) 
#     stop("min_mean must be greater than 0.")
#   if(max_mean_percentile < 0 || max_mean_percentile > 1)
#     stop("max_mean_percentile must be between 0 and 1.")
#   
#   means <- apply(reference_data, 2, mean)
#   max_mean <- quantile(means, max_mean_percentile)
#   reference_data <- reference_data[, (means >= min_mean) & (means <= max_mean)]
#   
#   index_zero <- NULL
#   index_ZI_genes <- NULL
#   if(!is.null(percent_ZI)) {
#     if(percent_ZI > 1 || percent_ZI < 0) {
#       stop("percent_ZI must be between 0 and 1.")
#     }
#     # index_ZI_genes <- which(apply(reference_data, 2, 
#     #                               function(x) mean(x == 0) > threshold_ZI))
#     index_ZI_genes <- which(apply(reference_data, 2, function(x) {
#       perc_zero <- mean(x == 0)
#       if(perc_zero < 0.05) {
#         return(FALSE)
#       } else if(mean(x) < 30) {
#         return(FALSE) 
#       } else if(var(x) == 0){
#         return(FALSE)
#       }
#       tryCatch(
#         {
#           m <- mean(x[x != 0]) # Used to calculate initial value.
#           fit_zinb <- fitdistrplus::fitdist(x, dzinb, 
#                                             start = list(size = m^2 / (var(x) - m), 
#                                                          mu = m, 
#                                                          rho = perc_zero / 2),
#                                             lower = c(1, 0, 0))
#           fit_nb <- fitdistrplus::fitdist(x, dzinb, 
#                                           start = list(size = m^2 / (var(x) - m), 
#                                                        mu = m),
#                                           lower = c(0, 0))
#           L <- dzinb(x, fit_zinb$estimate[1], fit_zinb$estimate[2], fit_zinb$estimate[3])
#           bic_zinb <- -2 * sum(log(L)) + log(length(x)) * 3 # Three parameters in model.
#           L <- dzinb(x, fit_nb$estimate[1], fit_nb$estimate[2])
#           bic_nb <- -2 * sum(log(L)) + log(length(x)) * 2 # Two parameters in model.
#           return(bic_zinb < bic_nb)
#         },
#         error = function(e) {
#           return(FALSE)
#         })
#     }))
#     index_ZI_genes <- unname(index_ZI_genes)
#     p_0 <- ceiling(p * percent_ZI)
#     if(p_0 > 0) {
#       if(length(index_ZI_genes) == 0) {
#         stop("None of the genes in the reference dataset are zero inflated.")
#       } else if(length(index_ZI_genes) < p_0) {
#         warning(paste("Not enough zero-inflated genes in reference dataset.",
#                       "Sampling with replacement will be used."))
#         index_zero <- sample(index_ZI_genes, p_0, replace = TRUE)
#       } else {
#         index_zero <- sample(index_ZI_genes, p_0)
#       }
#     }
#     p <- p - p_0
#   }
#   
#   index_nonzero <- NULL
#   if(p > 0) {
#     # If percent_zero_inflated is NULL, then index_genes will contain all genes.
#     # Otherwise, the zero inflated genes will be removed.
#     index_genes <- setdiff(1:ncol(reference_data), index_ZI_genes)
#     # Note, p will be zero if ceiling(p * percent_zero_inflated) is 1.
#     if(p > 0)
#       if(length(index_genes) == 0) {
#         stop("None of the genes in the reference dataset are non-zero-inflated.")
#       } else if(length(index_genes) < p) {
#         warning(paste("Not enough genes in reference dataset.",
#                       "Sampling with replacement will be used."))
#         index_nonzero <- sample(index_genes, p, replace = TRUE)
#       } else {
#         index_nonzero <- sample(index_genes, p)
#       }
#   }
#   
#   index <- c(index_zero, index_nonzero)
#   
#   return(reference_data[, index])
# }
