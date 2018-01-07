library(SimSeq)
library(SeqNet)
library(e1071)

validate_distribution_counts <- function(overdispersion = 10,
                                         intensity = 1.5,
                                         k = 1.5, 
                                         p = 1000, 
                                         seed) {
  
  cat("Seed:", seed, 
      "\t - overdispersion =", overdispersion,
      "\t intensity =", intensity, 
      "\t k =", k,
      "\n")
  set.seed(seed)
  data(kidney)
  ref <- t(kidney$counts) # Rows are samples, columns are genes.
  ref <- ref[, apply(ref, 2, mean) > 10] # Remove genes with mean < 30
  ref <- ref[, sample(ncol(ref))] # Mix the columns
  
  ref <- ref[, 1:p]
  
  # network <- create_network(p = p,
  #                           cliques = NULL,
  #                           hubs = NULL,
  #                           modules = NULL) # Network with no structure
  n <- nrow(ref)
  mu <- apply(ref, 2, mean)
  # x <- gen_gamma_poisson(n, network, mu, overdispersion, intensity, k)$x
  
  # par(mfrow = c(5, 2))
  # for(i in c(1, 11, 16, 91, 92)) {
  #   hist(x[, i], main = paste("x: mu =", round(mean(x[, i]), 2), "\n",
  #                             "kurt =", round(kurtosis(x[, i]), 2)),
  #        xlim = c(min(x[, i], ref[, i]), max(x[, i], ref[, i])))
  #   hist(ref[, i], main = paste("ref: mu =", round(mean(ref[, i]), 2), "\n",
  #                               "kurt =", round(kurtosis(ref[, i]), 2)),
  #        xlim = c(min(x[, i], ref[, i]), max(x[, i], ref[, i])))
  # }
  # par(mfrow = c(2, 2))
  # kurt <- apply(cbind(x, ref), 2, kurtosis)
  # avg <- apply(cbind(x, ref), 2, mean)
  # hist(apply(x, 2, kurtosis), 
  #      ylim = c(0, 0.7 * p),
  #      xlim = c(min(kurt), max(kurt)))
  # hist(apply(ref, 2, kurtosis), 
  #      ylim = c(0, 0.7 * p),
  #      xlim = c(min(kurt), max(kurt)))
  # plot(apply(x, 2, kurtosis), apply(x, 2, mean), 
  #      ylim = c(min(avg), max(avg)), xlim = c(min(kurt), max(kurt)))
  # plot(apply(ref, 2, kurtosis), apply(ref, 2, mean), 
  #      ylim = c(min(avg), max(avg)), xlim = c(min(kurt), max(kurt)))
  # 
  
  perc <- 0.20
  # par(mfrow = c(1, 1))
  network <- create_network(p = p,
                            cliques = NULL, #list(1:10, 11:20, 21:30, 31:40),
                            hubs = NULL, #list(41:50, 51:60, 61:70, 71:80, 81:100, 101:200),
                            modules = list(sample(1:p, p * perc), 
                                           sample(1:p, p * perc),
                                           sample(1:p, p * perc),
                                           sample(1:p, p * perc)), 
                            module_prob = 0.10)
  # plot(network, as_subgraph = TRUE)
  x <- gen_gamma_poisson(n, network, mu, overdispersion, intensity, k)$x
  
  # par(mfrow = c(5, 2))
  # for(i in c(1, 11, 16, 91, 92)) {
  #   hist(x[, i], main = paste("x: mu =", round(mean(x[, i]), 2), "\n",
  #                             "kurt =", round(kurtosis(x[, i]), 2)),
  #        xlim = c(min(x[, i], ref[, i]), max(x[, i], ref[, i])))
  #   hist(ref[, i], main = paste("ref: mu =", round(mean(ref[, i]), 2), "\n",
  #                               "kurt =", round(kurtosis(ref[, i]), 2)),
  #        xlim = c(min(x[, i], ref[, i]), max(x[, i], ref[, i])))
  # }
  # par(mfrow = c(2, 2))
  kurt <- apply(cbind(x, ref), 2, kurtosis)
  avg <- apply(cbind(x, ref), 2, mean)
  # hist(apply(x, 2, kurtosis), 
  #      ylim = c(0, 0.7 * p),
  #      xlim = c(min(kurt), max(kurt)))
  # hist(apply(ref, 2, kurtosis), 
  #      ylim = c(0, 0.7 * p),
  #      xlim = c(min(kurt), max(kurt)))
  if(overdispersion == param_set$overdispersion[1] &
     intensity == param_set$intensity[1] &
     k == param_set$k[1]) {
    plot(apply(ref, 2, kurtosis), apply(ref, 2, mean), 
         main = paste0("Reference\np = ", p, " (", perc * 100, "%)"),
         ylim = c(min(avg), max(avg)), xlim = c(min(kurt), max(kurt)))
  }
  plot(apply(x, 2, kurtosis), apply(x, 2, mean), 
       main = paste("overdis =", overdispersion, "\nintensity =", intensity, "\nk =", k),
       ylim = c(min(avg), max(avg)), xlim = c(min(kurt), max(kurt)))
}

param_set <- expand.grid(c(50, 60), 
                         c(1.5, 1.6),
                         c(1.5, 1.6))
colnames(param_set) <- c("overdispersion", "intensity", "k")
param_set <- rbind(c(100, 1.5, 1.5), 
                   c(70, 1.5, 1.5),
                   c(70, 1.6, 1.6), 
                   param_set)
p <- 5000
seed <- sample(1:10000, 1)

par(mfrow = c(3, 4))
apply(param_set, 1, function(param) {
      validate_distribution_counts(overdispersion = param["overdispersion"], 
                                   intensity = param["intensity"], 
                                   k = param["k"],
                                   p = p,
                                   seed = seed)
})
