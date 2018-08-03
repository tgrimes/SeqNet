#' @title Get performance measures for network inference
#' @description Computes the sensitivity, specificity, f1score, tdr, and tndr
#'   to measure the accuracy of an inferred network on the true network.
#' @param pred a matrix indicating the predicted connections
#'   between the genes. Can be adjacency matrix or sparse association matrix 
#'   (i.e. with non-significant associations removed); all non-zero entries are
#'   changed to 1. 
#' @param true a network object or an adjacency matrix indicating the 
#' true underlying connections in the graph.
#' @return A list containing sensitivity, specificity, f1_score, tdr (true 
#' discovery rate), and tndr (true non-discovery rate). If the predicted network
#' is empty, then tdr is set to NA. Similarly, if the predicted network is completely
#' connected, tndr is set to NA.
#' @export
get_network_inference_performance <- function(pred, true) {
  warning("depreciated: use get_confusion_matrix() instead.")
  
  if(class(true) == "network") {
    true <- get_adj_matrix_from_network(true)
  } else if((!is.matrix(true))) {
    stop("true should be a network or matrix.")
  }
  if(!is.matrix(pred)) {
    stop("pred is not a matrix.")
  }
  if(dim(pred)[1] != dim(pred)[2]) {
    stop("pred is not a square matrix.")
  }
  if(dim(true)[1] != dim(true)[2]) {
    stop("true is not a square matrix.")
  }
  if(dim(pred)[1] != dim(true)[1]) {
    stop("true and pred are not of equal dimension.")
  }
  if(!isSymmetric(unname(true))) {
    stop("true is not symmetric.")
  }
  if(!isSymmetric(unname(pred))) {
    stop("pred is not symmetric.")
  }
  
  # Change any non-zero entries in pred to 1. 
  pred[pred != 0] <- 1
  
  # Take values in lower triangle.
  true <- true[lower.tri(true)]
  pred <- pred[lower.tri(pred)]
  
  # Setup:
  #           pred
  #           1  0
  #           ----
  # true  1 | A  B
  #       0 | C  D
  #
  # sensitivity = A / (A + B)
  # specificity = D / (C + D)
  # f1-score = harmonic mean of sensitivity and specificity
  # tdr = A / (A + C)
  # tndr = D / (B + D)
  
  pred_1 <- pred == 1 
  true_1 <- true == 1
  
  A <- sum(pred_1 & true_1)
  B <- sum(!pred_1 & true_1)
  C <- sum(pred_1 & !true_1)
  D <- sum(!pred_1 & !true_1)
  
  sensitivity <- A / (A + B)
  specificity <- D / (C + D)
  f1_score <-  2 / (1 / sensitivity + 1 / specificity)
  tdr <- ifelse((A + C) == 0, NA, A / (A + C))
  tndr <- ifelse((B + D) == 0, NA, D / (B + D))
  
  return(list(sensitivity = sensitivity, specificity = specificity,
              f1_score = f1_score, tdr = tdr, tndr = tndr))
}


#' Get performance measures for network inference
#' 
#' Computes the confusion matrix for an inferred network based on the true 
#' network.
#' @param true a network object or an adjacency matrix indicating the 
#' true underlying connections in the graph.
#' @param pred a matrix indicating the predicted connections
#' between the genes. Can be adjacency matrix or sparse association matrix 
#' (i.e. with non-significant associations removed); all non-zero entries are
#' changed to 1. 
#' @return The confusion matrix.
#' @export
get_confusion_matrix <- function(true, pred) {
  if(class(true) == "network") {
    true <- get_adj_matrix_from_network(true)
  } else if((!is.matrix(true))) {
    stop("true should be a network or matrix.")
  }
  if(!is.matrix(pred)) {
    stop("pred is not a matrix.")
  }
  if(dim(pred)[1] != dim(pred)[2]) {
    stop("pred is not a square matrix.")
  }
  if(dim(true)[1] != dim(true)[2]) {
    stop("true is not a square matrix.")
  }
  if(dim(pred)[1] != dim(true)[1]) {
    stop("true and pred are not of equal dimension.")
  }
  if(!isSymmetric(unname(true))) {
    stop("true is not symmetric.")
  }
  if(!isSymmetric(unname(pred))) {
    stop("pred is not symmetric.")
  }
  
  # Change any non-zero entries in pred to 1. 
  pred[pred != 0] <- 1
  
  # Take values in lower triangle.
  true <- true[lower.tri(true)]
  pred <- pred[lower.tri(pred)]
  
  # Setup:
  #           pred
  #           1  0
  #           ----
  # true  1 | A  B
  #       0 | C  D
  #
  # sensitivity = A / (A + B)
  # specificity = D / (C + D)
  # f1-score = harmonic mean of sensitivity and specificity
  # tdr = A / (A + C)
  # tndr = D / (B + D)
  
  pred_1 <- pred == 1 
  true_1 <- true == 1
  
  A <- sum(pred_1 & true_1)
  B <- sum(!pred_1 & true_1)
  C <- sum(pred_1 & !true_1)
  D <- sum(!pred_1 & !true_1)
  
  confusion <- matrix(c(A, B, C, D), 2, 2, byrow = TRUE)
  colnames(confusion) <- c("1", "0")
  rownames(confusion) <- c("1", "0")
  dimnames(confusion) <- list(true = c("1", "0"), 
                              pred = c("1", "0"))
  
  return(confusion)
}

#' Provides summary measures for confusion matrix
#' 
#' @param confusion The confusion matrix.
#' @return A list containing the sensitivity (i.e. true positive rate; recall), 
#' specificity (i.e. true negative rate), true discovery rate (i.e. positive 
#' predictive value; precision), true non-discovery rate (i.e. negative 
#' predictive value), and f1-score (geometric mean of tdr and sensitivity).
#' @export
summarize_confusion <- function(confusion) {
  A <- confusion[1, 1]
  B <- confusion[1, 2]
  C <- confusion[2, 1]
  D <- confusion[2, 2]
  
  sensitivity <- A / (A + B) # True positive rate (TPR); Recall.
  specificity <- D / (C + D) # True negative rate (TNR).
  tdr <- ifelse((A + C) == 0, NA, A / (A + C)) # Positive predictive value (PPV); Precision.
  tndr <- ifelse((B + D) == 0, NA, D / (B + D)) # Negative predictive value (NPV).
  f1_score <-  2 / (1 / tdr + 1 / sensitivity)
  
  return(list(sensitivity = sensitivity,
              specificity = specificity,
              tdr = tdr,
              tndr = tndr,
              f1_score = f1_score))
}

#Performs tests to determine if get_performance_measures() is working properly.
test_network_measures <- function() {
  PASS <- TRUE
  assoc_true <- matrix(c(0, 1, 1, 0, 0,
                         1, 0, 1, 0, 1,
                         1, 1, 0, 0, 0,
                         0, 0, 0, 0, 1,
                         0, 1, 0, 1, 0), 5, 5)
  
  assoc_test1 <- assoc_true
  measures_test1 <- unlist(get_network_inference_performance(assoc_test1, assoc_true))
  if(any(measures_test1 != rep(1, 5))) {
    warning("performance measures did not pass test 1.")
    PASS <- FALSE
  }
  
  assoc_test2 <- !assoc_true
  measures_test2 <- unlist(get_network_inference_performance(assoc_test2, assoc_true))
  if(any(measures_test2 != rep(0, 5))) {
    warning("performance measures did not pass test 2.")
    PASS <- FALSE
  }
  
  assoc_test3 <- matrix(c(0, 1, 0, 0, 0,
                          1, 0, 0, 0, 1,
                          0, 0, 0, 1, 1,
                          0, 0, 1, 0, 1,
                          0, 1, 1, 1, 0), 5, 5)
  measures_test3 <- unlist(get_network_inference_performance(assoc_test3, assoc_true))
  if(any(measures_test3 != rep(3/5, 5))) {
    warning("performance measures did not pass test 3.")
    PASS <- FALSE
  }
  
  assoc_test4 <- matrix(c(0, 1, 0, 0, 0,
                          1, 0, 0, 0, 1,
                          0, 0, 0, 1, 1,
                          0, 0, 1, 0, 1,
                          0, 1, 1, 1, 0), 5, 5) * rnorm(25, 0, 1)
  assoc_test4 <- t(assoc_test4) + assoc_test4
  measures_test4 <- unlist(get_network_inference_performance(assoc_test4, assoc_true))
  if(any(measures_test4 != rep(3/5, 5))) {
    warning("performance measures did not pass test 4.")
    PASS <- FALSE
  }
  
  
  if(PASS) {
    print("All tests were passed.")
  }
}