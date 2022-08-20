test_that("linear.stderr.expectation approximately matches simulated results", {
  n <- 100000
  f <- 0.25
  genotypes <- sample(0:2, n, replace = TRUE, prob = c(
    f^2,
    2 * f * (1 - f),
    (1 - f)^2
  ))
  phenotypes <- rnorm(n)
  observed <- summary(lm(phenotypes ~ genotypes))$coeff[2, 2]
  expected <- linear.stderr.expectation(f, n)
  expect_true(abs(observed - expected) / observed < 0.01)
})
