
# Required:
#   edgeR - for calcNormFactors() and cpm().
#   WGCNA - for adjacency.fromSimilarity()
# library(edgeR)
# library(WGCNA)

#' Wrapper for WGCNA method
#' 
#' Conducts co-expression analysis using WGCNA method. Preprocessing is performed
#' on x: the trimmed mean of M-values (TMM) normalization, which accounts 
#' for between-sample biases, and counts per million (CPM) scaling, which 
#' addresses any difference in library sizes.
#' @param x The n by p matrix of counts.
#' @param threshold Cutoff for significant associations. If NULL, all scores
#' are returned. Otherwise, scores at or below this threshold are set to zero. 
#' @param method The method used to compute correlations. Should be either "pearson"
#'  or "spearman". The default is Spearman, which provides 
#'  the more conservative estimation of associations.
#' @return A list containing the p by p matrix of association scores, and the 
#' threshold used to determine significant associations.
#' @export
run_wgcna <- function(x, threshold = NULL, method = "pearson") {
  if(!(method %in% c("pearson", "spearman"))) {
    stop('method should be one of c("pearson", "spearman").')
  }
  
  #Normalize the data to account for 1) between-sample biases (TMM) and
  # 2) differeing library sizes (cpm).
  x <- t(t(x) * edgeR::calcNormFactors(x, method = "TMM"))
  x <- edgeR::cpm(x)
  
  # To prevent highly expressed outliers from dominating the between-gene
  # correlations, we perform a log(x + 1) transformation.
  x <- log(x + 1)
  
  # Use unsigned to treat positive and negative associations equally.
  scores <- WGCNA::adjacency(x, 
                             type = "unsigned", 
                             corFnc = "cor",
                             corOptions = list(use = "p", method = method)) 
  
  diag(scores) <- 0
  
  if(!is.null(threshold)) {
    scores[scores <= threshold] <- 0
  }
  
  colnames(scores) <- colnames(x)
  
  return(list(scores = scores, threshold = threshold))
}
