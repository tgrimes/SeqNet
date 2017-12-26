validate_distribution_counts <- function() {
  library(SimSeq)
  library(SeqNet)
  library(e1071)
  
  set.seed(1033)
  data(kidney)
  ref <- t(kidney$counts) # Rows are samples, columns are genes.
  ref <- ref[, apply(ref, 2, mean) > 10] # Remove genes with mean < 30
  ref <- ref[, sample(ncol(ref))] # Mix the columns
  rm(kidney)
  gc()
  
  p <- 1000
  ref <- ref[, 1:p]
  
  network <- create_network(p = p,
                            cliques = NULL,
                            hubs = NULL,
                            modules = NULL) # Network with no structure
  n <- nrow(ref)
  mu <- apply(ref, 2, mean)
  overdispersion <- 1
  intensity <- 1
  k <- 1.5
  x <- gen_gamma_poisson(n, network, mu, overdispersion, intensity, k)$x
  
  par(mfrow = c(5, 2))
  for(i in c(1, 11, 16, 91, 92)) {
    hist(x[, i], main = paste("x: mu =", round(mean(x[, i]), 2), "\n",
                              "kurt =", round(kurtosis(x[, i]), 2)),
         xlim = c(min(x[, i], ref[, i]), max(x[, i], ref[, i])))
    hist(ref[, i], main = paste("ref: mu =", round(mean(ref[, i]), 2), "\n",
                                "kurt =", round(kurtosis(ref[, i]), 2)),
         xlim = c(min(x[, i], ref[, i]), max(x[, i], ref[, i])))
  }
  par(mfrow = c(2, 2))
  kurt <- apply(cbind(x, ref), 2, kurtosis)
  avg <- apply(cbind(x, ref), 2, mean)
  hist(apply(x, 2, kurtosis), xlim = c(min(kurt), max(kurt)))
  hist(apply(ref, 2, kurtosis), xlim = c(min(kurt), max(kurt)))
  plot(apply(x, 2, kurtosis), apply(x, 2, mean), 
       ylim = c(min(avg), max(avg)), xlim = c(min(kurt), max(kurt)))
  plot(apply(ref, 2, kurtosis), apply(ref, 2, mean), 
       ylim = c(min(avg), max(avg)), xlim = c(min(kurt), max(kurt)))
  
  
  par(mfrow = c(1, 1))
  network <- create_network(p = p,
                            cliques = NULL,
                            hubs = NULL,
                            modules = list(1:200, 150:250, 240:400, 
                                           401:451, 451:501, 501:600), 
                            module_prob = 0.15)
  plot(network, as_subgraph = TRUE)
  x <- gen_gamma_poisson(n, network, mu, overdispersion, intensity, k)$x
  
  par(mfrow = c(5, 2))
  for(i in c(1, 11, 16, 91, 92)) {
    hist(x[, i], main = paste("x: mu =", round(mean(x[, i]), 2), "\n",
                              "kurt =", round(kurtosis(x[, i]), 2)),
         xlim = c(min(x[, i], ref[, i]), max(x[, i], ref[, i])))
    hist(ref[, i], main = paste("ref: mu =", round(mean(ref[, i]), 2), "\n",
                                "kurt =", round(kurtosis(ref[, i]), 2)),
         xlim = c(min(x[, i], ref[, i]), max(x[, i], ref[, i])))
  }
  par(mfrow = c(2, 2))
  kurt <- apply(cbind(x, ref), 2, kurtosis)
  avg <- apply(cbind(x, ref), 2, mean)
  hist(apply(x, 2, kurtosis), xlim = c(min(kurt), max(kurt)))
  hist(apply(ref, 2, kurtosis), xlim = c(min(kurt), max(kurt)))
  plot(apply(x, 2, kurtosis), apply(x, 2, mean), 
       ylim = c(min(avg), max(avg)), xlim = c(min(kurt), max(kurt)))
  plot(apply(ref, 2, kurtosis), apply(ref, 2, mean), 
       ylim = c(min(avg), max(avg)), xlim = c(min(kurt), max(kurt)))
}