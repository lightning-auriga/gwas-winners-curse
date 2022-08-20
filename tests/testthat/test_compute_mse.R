test_that("compute.beta.mse returns beta unadjusted when debiased and biased estimates are identical", {
  freq <- 0.25
  n <- 10000
  biased.beta <- 0.1
  ## for linear regression, se isn't really "biased" per se,
  ## but that's ok. value is approximately 0.01632993
  biased.se <- linear.stderr.expectation(freq, n)
  unbiased.beta <- biased.beta
  ## result should be unchanged
  observed <- compute.beta.mse(biased.beta, biased.se, unbiased.beta)
  expect_equal(observed, unbiased.beta)
})

## this test is one of many unfortunate tests in this package, as the true
## values are largely unknown. we're mostly judging this based on *consistency*
## between attempts and implementations, which is less than ideal.
test_that("compute.beta.mse returns appropriate mse adjustment when input betas differ", {
  freq <- 0.25
  n <- 10000
  biased.beta <- 0.1
  ## for linear regression, se isn't really "biased" per se,
  ## but that's ok. value is approximately 0.01632993
  biased.se <- linear.stderr.expectation(freq, n)
  ## correction does not require that the unbiased beta actually be the
  ## solution of debiasing
  unbiased.beta <- 0
  ## result should be different
  expected <- 0.002597402
  observed <- compute.beta.mse(biased.beta, biased.se, unbiased.beta)
  expect_equal(observed, expected, tolerance = 1e-5)
})

test_that("compute.ci.mse returns ci unadjusted when debiased and biased estimates are identical", {
  freq <- 0.25
  n <- 10000
  biased.beta <- 0.1
  ## value is approximately 0.01632993
  biased.se <- linear.stderr.expectation(freq, n)
  unbiased.beta <- biased.beta
  p.z <- 0.025
  unbiased.lci <- calculate.ci(unbiased.beta, biased.se, p.z)
  ## result should be unchanged
  observed <- compute.ci.mse(biased.beta, biased.se, unbiased.lci, p.z)
  expect_equal(observed, unbiased.lci)
})

## this test is one of many unfortunate tests in this package, as the true
## values are largely unknown. we're mostly judging this based on *consistency*
## between attempts and implementations, which is less than ideal.
test_that("compute.ci.mse returns appropriate mse adjustment when input cis differ", {
  freq <- 0.25
  n <- 10000
  biased.beta <- 0.1
  ## value is approximately 0.01632993
  biased.se <- linear.stderr.expectation(freq, n)
  unbiased.beta <- 0
  p.z <- 0.025
  ## -0.03200607
  unbiased.lci <- calculate.ci(unbiased.beta, biased.se, p.z)
  ## result should be different
  expected <- -0.02940867
  observed <- compute.ci.mse(biased.beta, biased.se, unbiased.lci, p.z)
  expect_equal(observed, expected, tolerance = 1e-5)
})
