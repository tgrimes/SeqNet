
# Required:
#   edgeR - for calcNormFactors() and cpm().
#   glasso - for adjacency.fromSimilarity()
# library(edgeR)
# library(huge)

#' Wrapper for glasso method
#' 
#' Conducts co-expression analysis using glasso method. P
#' @param x The n by p matrix of counts.
#' @param threshold Cutoff for significant associations. If NULL, all scores
#' are returned. Otherwise, scores at or below this threshold are set to zero. 
#' @param method The method used to compute correlations. Should be either "pearson"
#'  or "spearman". The default is Spearman, which provides 
#'  the more conservative estimation of associations.
#' @return A list containing the p by p matrix of association scores, and the 
#' threshold used to determine significant associations.
#' @export
run_glasso <- function(x, 
                       method = c("mb", "glasso", "ct"), 
                       criterion = c("ric", "stars"), 
                       verbose = FALSE, ...) {
  method <- tolower(method[1])
  criterion <- tolower(criterion[1])
  if(!(method %in% c("mb", "glasso", "ct"))) {
    stop('method should be one of "mb", "glasso", or "ct".')
  }
  if(criterion != "ric" & criterion != "stars") {
    stop("criterion must be either \"ric\" or \"stars\".")
  }
  
  # Use unsigned to treat positive and negative associations equally.
  x_huge <- huge(x, method = method, verbose = verbose)
  result <- huge.select(x_huge, criterion = criterion, verbose = verbose)
  
  # Choose largest lambda at or below the optimal value.
  index <- which(result$lambda <= result$opt.lambda)[1]
  if(length(index) == 0) index <- 10
  
  if(verbose) cat("Optimal lambda:", result$opt.lambda, "\n")
  
  scores <- as.matrix(x_huge$icov[[index]])
  diag(scores) <- 0
  
  colnames(scores) <- colnames(x)
  
  return(list(scores = scores))
}






# function(x, criterion = "ric") {
#   if(criterion != "ric" & criterion != "stars") {
#     stop("criterion must be either \"ric\" or \"stars\".")
#   }
#   
#   #Normalize the data to account for 1) between-sample biases (TMM) and
#   # 2) differeing library sizes (cpm).
#   x <- t(t(x) * calcNormFactors(x, method = "TMM"))
#   x <- cpm(x)
#   
#   #To prevent highly expressed outliers from dominating the between-gene
#   # correlations, we perform a log(x + 1) transformation.
#   x <- log(x + 1)
#   
#   x_huge <- huge(x, method = "glasso")
#   result <- huge.select(x_huge, criterion = criterion)
#   adj_matrix <- as.matrix(result$path[[1]])
#   
#   return(list(adj_matrix = adj_matrix))
# }
