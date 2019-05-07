testthat::test_that("qzinb works properly.", {
  n <- 100
  p <- 100
  params <- sapply(1:p, function(i) c(rnorm(1, 10, 2), rnorm(1, 500, 100), 0))
  vals <- apply(params, 2, function(x) rzinb(n, x[1], x[2], x[3]))
  x <- sapply(1:p, function(i) pzinb(vals[, i], params[1, i], 
                                     params[2, i], params[3, i]))
  temp1 <- x + 0.0
  temp2 <- x + 0.0
  
  overwrite_with_qzinb_cpp(temp1, params[1, ], params[2, ], params[3, ])
  for(i in 1:p) {
    temp2[, i] <- qzinb(temp2[, i], 
                        size = params[1, i], 
                        mu = params[2, i],
                        rho = params[3, i])
  }
  
  expect_true(all(temp1 == vals))
  expect_true(all(temp2 == vals))
})