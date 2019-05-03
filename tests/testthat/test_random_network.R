testthat::test_that("networks structures can be created", {
  network <- random_network(100)
  expect_is(network, "network")
  expect_is(network$modules[[1]], "network_module")
})