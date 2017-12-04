
# Required:
#   edgeR - for calcNormFactors() and cpm().
#   WGCNA - for adjacency.fromSimilarity()
# library(edgeR)
# library(WGCNA)

#' Wrapper for WGCNA method
#' 
#' Conducts co-expression analysis using WGCNA method.
#' @param x The n by p matrix of counts
#' @param threshold Cutoff for significant associations
#' @param beta Tuning parameter for WGCNA
#' @return A list containing `scores`, a p by p matrix of association scores, and 
#' `adj_matrix`, a p by p adjacency matrix.
#' @export
run_wgcna <- function(x, threshold = 0, beta = 6) {
  #Normalize the data to account for 1) between-sample biases (TMM) and
  # 2) differeing library sizes (cpm).
  x <- t(t(x) * edgeR::calcNormFactors(x, method = "TMM"))
  x <- edgeR::cpm(x)
  
  #To prevent highly expressed outliers from dominating the between-gene
  # correlations, we perform a log(x + 1) transformation.
  x <- log(x + 1)
  
  #Spearmans's correlation is used.
  x <- cor(x, method = "spearman")
  x[which(is.na(x))] <- 0
  
  assoc_matrix <- WGCNA::adjacency.fromSimilarity(x, power = beta)
  diag(assoc_matrix) <- 0
  adj_matrix <- 1 * (round(assoc_matrix, 2) > threshold)
  
  return(list(scores = assoc_matrix, adj_matrix = adj_matrix))
}