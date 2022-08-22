relative.tol <- 0.01

test_that("compute.variance.explained works based on linear simplification", {
  beta <- 0.1
  f <- 0.22
  expected <- 0.003432
  observed <- compute.variance.explained(beta, f)
  expect_equal(observed, expected)
})

test_that("compute.power works for simple linear model", {
  beta <- 0.1
  f <- 0.22
  n <- 10000
  p.thresh <- 1e-5
  expected <- 0.9252294
  observed <- compute.power(beta, f, n, p.thresh)
  expect_true((observed - expected) / expected < relative.tol)
})

test_that("expected.replication.count sums power across multiple variants", {
  beta <- c(0.1, 0.2)
  f <- c(0.22, 0.12)
  n <- c(10000, 5000)
  p.thresh <- c(1e-5, 1e-6)
  expected <- 1.871267
  observed <- expected.replication.count(beta, f, n, p.thresh)
  expect_true((observed - expected) / expected < relative.tol)
})
