#' Sample genes from reference dataset
#' 
#' @param reference_data The reference data.frame to use.
#' @param p The number of genes (columns) to sample.
#' @param percent_ZI The desired percentage of zero-inflated genes. This 
#' percentage of zero-inflated genes will be sampled from the reference
#' dataset, and the remaining will be non-zero-inflated. If
#' NULL, then genes are sampled at random from the reference dataset.
#' @param threshold_ZI The minimum proportion of zero counts for a gene to be
#' considered as zero inflated. This is used to identify which genes in the
#' reference dataset are zero-inflated.
#' @return The modified reference dataset.
#' @note If p is greater than the number of columns in the reference dataset,
#' then sampling with replacement will be used (with a warning message). 
#' @export 
#' @examples 
#' \donttest{
#' data(reference)
#' rnaseq <- reference$rnaseq
#' rnaseq_subset <- sample_reference_data(rnaseq, 10)
#' }
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
                      "Sampling with replacement is used."))
        index_nonzero <- sample(index_genes, p, replace = TRUE)
      } else if(length(index_genes) > p) {
        index_nonzero <- sample(index_genes, p)
      } else {
        index_nonzero <- index_genes
      }
  }
  
  index <- c(index_zero, index_nonzero)
  
  return(reference_data[, index])
}