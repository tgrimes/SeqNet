
# Required:
#   edgeR - for calcNormFactors() and cpm().
#   WGCNA - for adjacency.fromSimilarity()
# library(edgeR)
# library(WGCNA)

#' Wrapper for WGCNA method
#' 
#' Conducts co-expression analysis using WGCNA method. 
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
