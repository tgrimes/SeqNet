testthat::test_that("counts can be generated", {
  reference <- matrix(rpois(200, 10), nrow = 10, ncol = 20)
  x <- gen_rnaseq(10, gen_partial_correlations(create_empty_network(20)), reference)$x
  expect_is(x, "matrix")
  x <- gen_rnaseq(10, gen_partial_correlations(random_network(20)), reference)$x
  expect_is(x, "matrix")
  x <- gen_rnaseq(20, gen_partial_correlations(random_network(20)), reference)$x
  expect_is(x, "matrix")
  x <- gen_rnaseq(20, gen_partial_correlations(random_network(10)), reference)$x
  expect_is(x, "matrix")
})