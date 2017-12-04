# Example code and testing.

#' Example of SeqNet package
#' 
#' @export
run_example <- function() {
  set.seed(6720)
  p <- 100
  
  network <- create_network(p = p,
                            clique_size = c(5, 5),
                            hub_size = c(5, 5, 4, 6, 6, 6),
                            module_size = c(35),
                            module_prob = 0.03,
                            nonoverlapping = TRUE)
  plot(network)

  network <- create_network(p = p,
                            clique_size = c(5, 5),
                            hub_size = c(5, 5, 4, 6, 6, 6),
                            module_size = c(35),
                            module_prob = 0.03,
                            nonoverlapping = FALSE)
  plot(network)
  
  network <- create_network(p = p,
                            cliques = list(1:5, 6:10),
                            hubs = list(11:15, 16:25, 
                                        26:29, 
                                        c(27, 31:35), 
                                        c(28, 36:40),
                                        c(29, 41:45)),
                            modules = list(56:90), 
                            module_prob = 0.03)
  plot(network)
  
  network <- add_sign_to_network(network)
  mu <- sample(get_reference_count_means(), p, replace = TRUE)
  
  n <- 50
  overdispersion_param <- 1
  intensity_param <- 1
  k <- 1.5
  df <- gen_gamma_poisson(n, network, mu, 
                          overdispersion_param = overdispersion_param, 
                          intensity_param = intensity_param,
                          k = k)
  # Contents of the list returned from the generator:
  #   df$x - the generated sample as an n by p matrix. 
  #   df$network - the underlying regulatory network used to generate the sample.
  #   df$mu - the default average expression value for each gene.
  
  cpls <- run_cpls(df$x, v = 3)$scores
  corr <- run_corr(df$x)$scores
  wgcna <- run_wgcna(df$x)$scores
  
  n_top <- 100 # sum(get_adj_matrix_from_network(network))
  adj_matrix_using_top_scores <- function(scores) {
    scores[abs(scores) < sort(abs(scores), decreasing = TRUE)[n_top]] <- 0
    scores[scores != 0] <- 1
    return(scores)
  }
  
  par(mfrow = c(2, 2))
  g <- plot(network, main = "Regulatory Network")
  plot.network(adj_matrix_using_top_scores(cpls), g, main = "cPLS")
  plot.network(adj_matrix_using_top_scores(corr), g, main = "cor")
  plot.network(adj_matrix_using_top_scores(wgcna), g, main = "WGCNA")
  
  # View(df$x)
}