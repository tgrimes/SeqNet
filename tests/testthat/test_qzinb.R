testthat::test_that("qzinb works properly.", {
  n <- 100
  p <- 100
  params <- sapply(1:p, function(i) c(rnorm(1, 10, 2), rnorm(1, 500, 100), 0))
  vals <- apply(params, 2, function(x) rzinb(n, x[1], x[2], x[3]))
  x <- sapply(1:p, function(i) pzinb(vals[, i], params[1, i], 
                                     params[2, i], params[3, i]))
  vals_from_qzinb <- x + 0.0
  
  for(i in 1:p) {
    vals_from_qzinb[, i] <- qzinb(x[, i], 
                                  size = params[1, i], 
                                  mu = params[2, i],
                                  rho = params[3, i])
  }
  
  expect_true(all(vals_from_qzinb == vals))
})