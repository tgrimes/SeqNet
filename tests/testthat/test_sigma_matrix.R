testthat::test_that("sigma matrix is correct", {
  adj <- matrix(0, 5, 5)
  adj[1, -c(1)] <- 0.4
  adj[-c(1), 1] <- 0.4
  adj2 <- matrix(0, 3, 3)
  adj2[1, 2] <- 0.9
  adj2[2, 1] <- 0.9
  colnames(adj2) <- 6:8
  mod1 <- create_module_from_association_matrix(adj)
  mod2 <- create_module_from_association_matrix(adj2)
  nw <- create_network_from_modules(8, list(mod1, mod2))
  expect_true(nrow(mod2$edges) == 1)
  expect_true(ncol(mod2$edges) == 3)
  expect_is(nw, "network")
  expect_is(get_sigma(nw), "matrix")
  expect_true(all(get_sigma(nw)[8, -8] == 0))
})
