testthat::test_that("cpp functions work", {
  x <- rnorm(100)
  expect_equal(ecdf(x)(x), ecdf_cpp(x))
  x <- rep(1, 10)
  expect_equal(ecdf(x)(x), ecdf_cpp(x))
  x <- rep(10, 1)
  expect_equal(ecdf(x)(x), ecdf_cpp(x))
  x <- 1
  expect_equal(ecdf(x)(x), ecdf_cpp(x))
  
  adj <- matrix(0, 5, 5)
  edges <- rbind(c(1, 2), c(3, 4), c(3, 5))
  adj[edges] <- 1
  adj[edges[, 2:1]] <- 1
  expect_equal(components_in_adjacency(adj), c(1, 1, 2, 2, 2))
  
  expect_true(is_adjacency_cpp(adj))
  
  expect_equal(edges_from_adjacency_cpp(adj), edges)
  
  expect_equal(ring_lattice_cpp(4, 0), 
               matrix(0, 4, 4))
  expect_equal(ring_lattice_cpp(4, 1), 
               matrix(c(0, 1, 0, 1, 1, 0, 1, 0,
                        0, 1, 0, 1, 1, 0, 1, 0), 4, 4))
  mat <- matrix(1, 4, 4)
  diag(mat) <- 0
  expect_equal(ring_lattice_cpp(4, 2), mat)
})