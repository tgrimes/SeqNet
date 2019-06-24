testthat::test_that("networks can be created from matricies", {
  adj <- random_module_structure(10)
  colnames(adj) <- 1:10
  expect_is(adj, "matrix")
  expect_true(is_adjacency_cpp(adj))
  network <- create_network_from_adjacency_matrix(adj)
  expect_is(network, "network")
  expect_equal(get_adjacency_matrix(network), adj)
  
  network <- gen_partial_correlations(network)
  assoc <- get_association_matrix(network)
  expect_is(assoc, "matrix")
  network <- create_network_from_association_matrix(assoc)
  expect_is(network, "network")
  expect_equal(get_association_matrix(network), assoc)
  expect_equal(get_adjacency_matrix(network), get_adjacency_matrix(assoc))
})