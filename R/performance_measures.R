#' Get confusion matrix for predicted classes.
#' 
#' Computes the confusion matrix for predicted classes. The current implementation
#' only handles binary classes and is designed to primarily deal with network objects
#' and adjacency matricies. The true values are assumed to be either 0 or 1,
#' and the predicted values are 0 or nonzero - all nonzero predicted values are 
#' converted to 1.
#' @param true a network object, adjacency matrix indicating the 
#' true underlying connections in the graph, or a binary vector.
#' @param pred a matrix indicating the predicted connections
#' between the genes. Can be adjacency matrix or sparse association matrix 
#' (i.e. with non-significant associations set to zero); all non-zero entries 
#' will be set to 1. If true is a binary vector, pred can also be a vector.
#' @return An object of class "CONFUSION".
#' @export
get_confusion_matrix <- function(true, pred) {
  if(class(true) == "network") {
    true <- get_adj_matrix_from_network(true)
    if(!is.matrix(pred)) {
      stop("true is a network object, so pred should be a matrix.")
    }
  } else if(is.matrix(true) & !is.matrix(pred)) {
    stop("true is a matrix, but pred is not.")
  } else if(!is.matrix(true) & is.matrix(pred)) {
    stop("pred is a matrix, but true is not.")
  } 
  if(!all(true %in% c(0, 1)))
    stop("true should contain only 0's and 1's")
  
  # Change any non-zero entries in pred to 1. 
  pred[pred != 0] <- 1
  
  if(is.matrix(true) & is.matrix(pred)) {
    # Take values in lower triangle.
    true <- true[lower.tri(true)]
    pred <- pred[lower.tri(pred)]
  }
  
  if(length(true) != length(pred)) {
    stop("length of true is not equal to length of true")
  }
  
  # Setup:
  #             True value
  #                 1   0
  #                 ------
  # Predicted   1 | TP  FP
  #     Value   0 | FN  TN
  #
  # sensitivity = TP / (TP + FN)
  # specificity = TN / (FP + TN)
  # tdr = TP / (TP + FP)
  # tndr = TN / (TN + FN)
  
  pred_1 <- pred == 1 
  true_1 <- true == 1
  
  TP <- sum(pred_1 & true_1)
  FN <- sum(!pred_1 & true_1)
  FP <- sum(pred_1 & !true_1)
  TN <- sum(!pred_1 & !true_1)
  
  mat <- matrix(c(TP, FN, FP, TN), 2, 2)
  rownames(mat) <- c("1", "0")
  colnames(mat) <- c("1", "0")
  dimnames(mat) <- list(`Predicted value` = c("1", "0"), 
                              `True value` = c("1", "0"))
  
  confusion <- list(mat = mat,
                    TP = TP,
                    FN = FN,
                    FP = FP,
                    TN = TN)
  class(confusion) <- "CONFUSION"
  return(confusion)
}

#' Coerce a matrix into "CONFUSION" object
#' 
#' Coerces a 2 by 2 matrix containing integer counts into an object of class
#' "CONFUSION".
#' @param mat A 2 by 2 matrix containing integer values.
#' @return An object of class "CONFUSION" containing the same counts as mat.
#' @export
as.CONFUSION <- function(mat) {
  if(!is.matrix(mat)) 
    stop("mat should be a matrix.")
  if(length(dim(mat)) != 2) 
    stop("mat does not have two dimensions.")
  if(!all(dim(mat) == 2)) 
    stop("mat is not a 2 by 2 matrix.")
  if(!all(sapply(as.vector(mat), is.integer))) 
    stop("mat must contain integer values.")
  
  rownames(mat) <- c("1", "0")
  colnames(mat) <- c("1", "0")
  dimnames(mat) <- list(`Predicted value` = c("1", "0"), 
                        `True value` = c("1", "0"))
  
  confusion <- list(mat = mat,
                    TP = mat[1, 1],
                    FN = mat[2, 1],
                    FP = mat[1, 2],
                    TN = mat[2, 2])
  class(confusion) <- "CONFUSION"
  return(confusion)
}



#' Print function for objects of class "CONFUSION"
#' 
#' Print the confusion matrix.
#' @param confusion An object of class "CONFUSION".
#' @param ... Additional arguments to be passed into the print function.
#' @export
print.CONFUSION <- function(confusion, ...) {
  print(confusion$mat, ...)
}



#' Coerce "CONFUSION" object into a data.frame
#' 
#' Function to coerce an object of class "CONFUSION" into a data.frame.
#' @param confusion An object of class "CONFUSION".
#' @param ... Additional arguments; these are currently ignored.
#' @return Returns a data.frame object.
#' @export
as.data.frame.CONFUSION <- function(confusion, ...) {
  return(cbind(expand.grid(dimnames(confusion$mat)),
         counts = as.vector(confusion$mat)))
}

#' Plots confusion matrix
#' 
#' Plot method for objects of class "CONFUSION".
#' @param confusion An object of class "CONFUSION".
#' @param main Title for the plot.
#' @param ... Other graphics parameters; these are currently ignored.
#' @return The generated ggplot object is returned. This can be used for further 
#' editing.
#' @note The gradient used to color each cell is proportional to the sum of
#' adjacent cells with the current cell; for example, the gradient to color 
#' the TP cell is TP / (TP + FN + FP). This re-weighting helps reveal good 
#' performance in the case of imbalanced class labels.
#' @export
plot.CONFUSION <- function(confusion, main = NULL, ...) {
  df <- as.data.frame(confusion)
  mat <- confusion$mat
  total <- sum(mat)
  w_counts <- c(mat[1, 1] / (total - mat[-1, -1]),
                mat[2, 1] / (total - mat[-2, -1]),
                mat[1, 2] / (total - mat[-1, -2]),
                mat[2, 2] / (total - mat[-2, -2]))
  labels <- c("TP", "FN", "FP", "TN")
  cbind(df, w_counts = w_counts, labels = labels)
  g <- ggplot(data = df, 
              aes(x = `True value`, y = `Predicted value`)) +
    scale_y_discrete(limits = rev(levels(df$`Predicted value`))) +
    scale_x_discrete(position = "top") +
    geom_tile(aes(fill = w_counts), color = "black") +
    geom_text(aes(label = paste(labels, "=", counts)), color = "white") +
    scale_fill_gradient(limits = c(0, max(sum(w_counts) / 2, max(w_counts)))) +
    theme_minimal() +
    theme(legend.position = "none")
  
  if(!is.null(main))
    g <- g + labs(title = main)
  
  print(g)
  return(g)
}


#' Provides summary measures for confusion matrix
#' 
#' @param confusion The confusion matrix. This can either be an object of class
#' "CONFUSION" or a 2 by 2 integer matrix.
#' @return A list containing several summary measures: 
#' sensitivity (a.k.a. true positive rate; recall), 
#' specificity (a.k.a. true negative rate), 
#' TDR (true discovery rate; a.k.a. positive predictive value; precision), 
#' TNDR (true non-discovery rate; a.k.a. negative predictive value), 
#' f1-score (harmonic mean of TDR and sensitivity), and
#' MCC (Matthew's correlation coefficient).
#' @export
summarize_confusion <- function(confusion) {
  if(class(confusion) != "CONFUSION") {
    confusion <- tryCatch(
      as.CONFUSION(confusion),
      error = function(e) {
        stop(paste("no method to hand object of class", class(confusion)))
      }
    )
  }
  
  TP <- confusion$TP
  FN <- confusion$FN
  FP <- confusion$FP
  TN <- confusion$TN
  
  sensitivity <- TP / (TP + FN) # True positive rate (TPR); Recall.
  specificity <- TN / (TN + FP) # True negative rate (TNR).
  TDR <- ifelse((TP + FP) == 0, NA, TP / (TP + FP)) # Positive predictive value (PPV); Precision.
  TNDR <- ifelse((FN + TN) == 0, NA, TN / (TN + FN)) # Negative predictive value (NPV).
  f1_score <-  2 / (1 / TDR + 1 / sensitivity)
  MCC <- (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  return(list(sensitivity = sensitivity,
              specificity = specificity,
              TDR = TDR,
              TNDR = TNDR,
              f1_score = f1_score,
              MCC = MCC))
}

#' Compute ROC curve 
#' 
#' @param true A vector of true labels.
#' @param pred A vector of predicted 
#' @return An object of class "ROC".
get_roc <- function(true, pred) {
  if(is.factor(true) & length(levels(true)) != 2)
    strop("true should only contain two factor levels.")
  if(!all(true %in% c(0, 1)))
    stop("true is not a factor nor is a vector of 0's and 1's.")
  
  if(is.matrix(pred)) pred <- pred[lower.tri(pred)]
  if(is.matrix(true)) true <- true[lower.tri(true)]
  
  ranked <- order(abs(pred), decreasing = TRUE)
  pred <- pred[ranked]
  true <- true[ranked]
  
  sens <- c(0, cumsum(true) / sum(true))
  spec <- c(1, 1 - cumsum(!true) / sum(!true))
  tdr <- c(0, cumsum(true) / (1:length(true)))
  tndr <- c(rev(cumsum(!true)) / (length(true):1), 0)
  df <- data.frame(sens = sens, 
                   spec = spec, 
                   tdr = tdr,
                   tndr = tndr)
  class(df) <- "ROC"
  return(df)
}

print.ROC <- function(roc, digits = 4, ...) {
  n <- length(roc$sens) - 1 
  auc <- get_auc(roc)
  max_f1_score <- 2 * prod(c(roc$tdr[n + 1], roc$sens[n + 1])) / 
    sum(c(roc$tdr[n + 1], roc$sens[n + 1]))
  tdr_range <- range(roc$tdr[-1])
  tndr_range <- range(roc$tndr[1:n])
  
  cat("Sample size:", n,
      "\nTrue discovery rate: from", round(tdr_range[1], digits), 
        "to", round(tdr_range[2], digits), 
      "\nTrue non-discovery rate: from", round(tndr_range[1], digits), 
        "to", round(tndr_range[2], digits), 
      "\nArea under Curve:", round(auc, digits),
      "\nMaximum F1 Score:", round(max_f1_score, digits), "\n")
}

plot_roc <- function(roc, title = "ROC curve", ...) {
  if(class(roc) != "ROC")
    stop(paste("not applicable for an object of class", class(roc)))
  
  auc <- get_auc(x = (1 - roc$spec), y = roc$sens)
  g <- data.frame(spec = roc$spec, sens = roc$sens) %>%
    as.data.frame() %>%
    ggplot(aes(x = 1 - spec, y = sens)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = 3) +
    theme_bw() +
    lims(y = c(0, 1),
         x = c(0, 1)) +
    labs(x = "1 - specificity", 
         y = "sensitivity",
         title = title) +
    geom_text(x = 0.8, y = 0.1, size = 6,
              aes(label = paste("AUC =", round(auc, 3))))
  plot(g)
  return(g)
}


plot_pr <- function(roc, title = "Precision-Recall curve", ...) {
  if(class(roc) != "ROC")
    stop(paste("not applicable for an object of class", class(roc)))
  
  auc <- get_auc(x = roc$sens, y = roc$tdr)
  
  g <- data.frame(sens = roc$sens, tdr = roc$tdr) %>%
    ggplot(aes(sens, tdr)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = 3) +
    theme_bw() +
    lims(y = c(0, 1),
         x = c(0, 1)) +
    labs(x = "Recall", 
         y = "Precision",
         title = title) +
    geom_text(x = 0.8, y = 0.1, size = 6,
              aes(label = paste("AUC =", round(auc, 3))))
  plot(g)
  return(g)
}

plot.ROC <- plot_roc

# Example for get_auc():
# get_auc(get_roc(rbinom(100, 1, 0.5), runif(100, 0, 1)))

#' Get area under the curve
#' 
#' Can be applied to an roc object to obtain AUROC, or the user may specify
#' values for x and y.
#' @param roc An object of class roc.
#' @param x An optional vector of x-coordinates for the curve.
#' @param y An optional vector of y-coordinates for the curve.
#' @return The estimated area under the curve.
get_auc <- function(roc = NULL, x = NULL, y = NULL) {
  if(!is.null(roc)) {
    if(class(roc) != "ROC")
      stop(paste("not applicable for an object of class", class(roc)))
    if(is.null(x))
      x <- 1 - roc$spec
    if(is.null(y))
      y <- roc$sens
  }
  
  if(length(x) != length(y)) 
    stop("x and y have differening lengths.")
  
  y <- y[-1]
  n <- length(x)
  auc <- sum(y * (x[-1] - x[-n]))
  return(auc)
}


# Performs tests to determine if get_confusion_matrix() and summarize_confusion()
# are working properly.
test_network_measures <- function() {
  PASS <- TRUE
  assoc_true <- matrix(c(0, 1, 1, 0, 0,
                         1, 0, 1, 0, 1,
                         1, 1, 0, 0, 0,
                         0, 0, 0, 0, 1,
                         0, 1, 0, 1, 0), 5, 5)
  
  assoc_test1 <- assoc_true
  
  measures_test1 <- unlist(
    summarize_confusion(get_confusion_matrix(assoc_true, assoc_test1)))
  if(any(measures_test1 != rep(1, 6))) {
    warning("performance measures did not pass test 1.")
    PASS <- FALSE
  }
  
  assoc_test2 <- !assoc_true
  measures_test2 <- unlist(
    summarize_confusion(get_confusion_matrix(assoc_true, assoc_test2)))  
  if(any(measures_test2[-6] != rep(0, 5)) || measures_test2[6] != -1) {
    warning("performance measures did not pass test 2.")
    PASS <- FALSE
  }
  
  assoc_test3 <- matrix(c(0, 1, 0, 0, 0,
                          1, 0, 0, 0, 1,
                          0, 0, 0, 1, 1,
                          0, 0, 1, 0, 1,
                          0, 1, 1, 1, 0), 5, 5)
  measures_test3 <- unlist(
    summarize_confusion(get_confusion_matrix(assoc_true, assoc_test3))) 
  if(any(measures_test3[-6] != rep(3/5, 5))) {
    warning("performance measures did not pass test 3.")
    PASS <- FALSE
  }
  
  assoc_test4 <- matrix(c(0, 1, 0, 0, 0,
                          1, 0, 0, 0, 1,
                          0, 0, 0, 1, 1,
                          0, 0, 1, 0, 1,
                          0, 1, 1, 1, 0), 5, 5) * rnorm(25, 0, 1)
  assoc_test4 <- t(assoc_test4) + assoc_test4
  measures_test4 <- unlist(
    summarize_confusion(get_confusion_matrix(assoc_true, assoc_test4))) 
  if(any(measures_test4[-6] != rep(3/5, 5))) {
    warning("performance measures did not pass test 4.")
    PASS <- FALSE
  }
  
  
  if(PASS) {
    print("All tests were passed.")
  }
}
