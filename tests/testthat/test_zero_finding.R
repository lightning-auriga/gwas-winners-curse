test_that("calculate.ci returns a valid confidence interval estimate assuming normality", {
  expected <- 126.6736
  observed <- calculate.ci(127.3, 0.31961, 0.025)
  expect_equal(observed, expected, tolerance = 1e-5)
})

test_that("calculate.adjusted.trait corrects a phenotype based on variant summary metrics", {
  trait.mean <- -0.5
  freq <- 0.25
  beta <- 0.1
  expected <- -0.55
  observed <- calculate.adjusted.trait(trait.mean, freq, beta)
  expect_equal(observed, expected, tolerance = 1e-5)
})

test_that("debiasing.func iterates a single step successfully", {
  beta.debiased <- 0
  p.thresh <- 1e-5
  beta.biased <- 0.1
  se.biased <- 0.01632993
  expected <- -0.1
  observed <- debiasing.func(beta.debiased, p.thresh, beta.biased, se.biased)
  expect_equal(observed, expected, tolerance = 1e-5)
})

test_that("debias.beta can debias a beta", {
  beta.biased <- -0.1191
  se.biased <- 0.0263
  freq <- 0.9092
  trait.mean <- 0
  p.thresh <- 1e-5
  expected <- -0.0123941
  observed <- debias.beta(beta.biased, se.biased, freq, trait.mean, p.thresh)
  expect_true(abs(observed - expected) / expected < 0.01)
})
