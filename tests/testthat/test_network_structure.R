testthat::test_that("networks have appropriate structure.", {
  network <- create_empty_network(10)
  adj <- get_adjacency_matrix(network)
  expect_true(all(adj == 0))
  
  network <- random_network(10, n_modules = 0)
  adj <- get_adjacency_matrix(network)
  expect_true(all(adj == 0))
  
  network <- random_network(10, neig_size = 1, min_module_size = 10, 
                            prob_rewire = 0, prob_remove = 0)
  vals <- get_network_characteristics(network)
  expect_equal(vals$p, 10)
  expect_equal(vals$n_edges, 10)
  expect_equal(vals$`avg degree`, 2)
  expect_equal(vals$`clustering coef`, 0)
  expect_equal(vals$`avg path length`, 25 / 9)
  
  network <- random_network(50, neig_size = 3, avg_module_size = 100,
                            min_module_size = 50, 
                            prob_rewire = 0, prob_remove = 0)
  vals <- get_network_characteristics(network)
  expect_equal(vals$p, 50)
  expect_equal(vals$n_edges, 150)
  expect_equal(vals$`avg degree`, 6)
  expect_equal(vals$`clustering coef`, 0.6)
  expect_equal(round(vals$`avg path length`, 3), 4.592)
})