test_that("poisson.binom.test functions as expected, greater mode", {
  expected <- pbinom(1, 2, 0.5, lower.tail = FALSE) * pbinom(2, 3, 0.45, lower.tail = FALSE)
  observed <- poisson.binom.test(5, c(0.5, 0.5, 0.45, 0.45, 0.45), alt = "greater")
  expect_equal(observed, expected)
})

test_that("poisson.binom.test functions as expected, less mode", {
  expected <- pbinom(0, 2, 0.5) * pbinom(0, 3, 0.45)
  observed <- poisson.binom.test(0, c(0.5, 0.5, 0.45, 0.45, 0.45), alt = "less")
  expect_equal(observed, expected)
})

test_that("poisson.binom.test functions as expected, two mode, at expectation", {
  p <- c(0.5, 0.5, 0.4, 0.4, 0.4, 0.4, 0.4)
  x <- 3
  expected <- 1
  observed <- poisson.binom.test(x, p, alt = "two")
  expect_equal(observed, expected)
})

test_that("poisson.binom.test functions as expected, two mode, below expectation", {
  expected <- pbinom(0, 2, 0.5) * pbinom(0, 3, 0.45) +
    pbinom(1, 2, 0.5, lower.tail = FALSE) * pbinom(2, 3, 0.45, lower.tail = FALSE)
  observed <- poisson.binom.test(0, c(0.5, 0.5, 0.45, 0.45, 0.45), alt = "two")
  expect_equal(observed, expected)
})

test_that("poisson.binom.test functions as expected, two mode, above expectation", {
  expected <- pbinom(1, 2, 0.5, lower.tail = FALSE) * pbinom(2, 3, 0.45, lower.tail = FALSE)
  observed <- poisson.binom.test(5, c(0.5, 0.5, 0.45, 0.45, 0.45), alt = "two")
  expect_equal(observed, expected)
})
