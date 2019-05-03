testthat::test_that("network module names are indexing correctly.", {
  network <- create_empty_network(20)
  network <- set_node_names(network, 1:20)
  module <- create_empty_module(11:20)
  struct <- random_module_structure(10)
  module <- set_module_edges(module, struct)
  network <- add_modules_to_network(network, module)
  adj <- get_adjacency_matrix(network)
  empty_struct <- matrix(0, 10, 10)
  colnames(empty_struct) <- 1:10
  module_struct <- get_adjacency_matrix(module)
  
  expect_equal(adj[1:10, 1:10], empty_struct)
  expect_equal(adj[11:20, 11:20], module_struct)
  
  network <- random_network(100, avg_module_size = 20, sd_module_size = 10,
                            n_modules = 5, consistent_connections = TRUE)
  adj <- get_adjacency_matrix(network)
  mod5 <- network$modules[[5]]
  adj_mod5 <- get_adjacency_matrix(mod5)
  index <- mod5$nodes
  expect_equal(adj[index, index], adj_mod5)
})