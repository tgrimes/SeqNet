
# Required:
#   edgeR - for calcNormFactors() and cpm().
#   WGCNA - for adjacency.fromSimilarity()
# library(edgeR)
# library(WGCNA)

#' Wrapper for WGCNA method
#' 
#' Conducts co-expression analysis using WGCNA method.
#' @param x The n by p matrix of counts.
#' @param threshold Cutoff for significant associations. Scores are rounded to
#' two decimal places. If threshold is provided, rounded scores at or below this 
#' threshold are set to zero. Otherwise, all scores are returned.
#' @param beta Tuning parameter for WGCNA.
#' @return A list containing the p by p matrix of association scores, and the 
#' threshold used to determine significant associations.
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
  
  scores <- WGCNA::adjacency.fromSimilarity(x, power = beta)
  diag(scores) <- 0
  
  if(!is.null(threshold)) {
    scores[round(scores, 2) <= threshold] <- 0
  }
  
  return(list(scores = scores, threshold = threshold))
}