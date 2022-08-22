#' @title
#' Performs an approximate test of a simple null hypothesis
#' about the probability of success in a Poisson binomial
#' series.
#'
#' @description
#' This function extends stats::binom.test for the Poisson
#' binomial distribution. Most of the original binom.test
#' functionality is removed in favor of a streamlined
#' calculation of the p-value itself.
#'
#' @details
#' To be entirely frank, I don't have a great deal of confidence
#' in whether this thing works or not. It seems like a clever trick
#' to me, but then, it's entirely possible it isn't implemented
#' in the source package itself for a reason.
#'
#' @param x Numeric number of success counts.
#' @param n Numeric number of trials.
#' @param p Numeric sequence of probabilities for each trial in turn.
#' @param alternative Character, one of "two.sided", "less", or "greater."
#' Analogous to the same parameter in binom.test.
#' @return Numeric probability of the alternative hypothesis
#' given the input parameters.
#' @export
#' @examples
#' x <- 20
#' n <- 30
#' p <- seq(0, 1, length.out = 10)
#' poisson.binom.test(x, n, p, "greater")
poisson.binom.test <- function(x, p = 0.5, alternative = c("two.sided", "less", "greater")) {
  stopifnot(length(x) == 1)
  stopifnot(is.numeric(x) || is.integer(x))
  stopifnot(x >= 0)
  alternative <- match.arg(alternative)
  p.value <- switch(alternative,
    less = PoissonBinomial::ppbinom(x, p),
    greater = PoissonBinomial::ppbinom(x - 1, p, lower.tail = FALSE),
    two.sided = {
      rel.err <- 1 + 1e-07
      d <- dpbinom(x, p)
      m <- sum(p)
      if (x == m) {
        1
      } else if (x < m) {
        i <- seq.int(from = ceiling(m), to = length(p))
        y <- sum(PoissonBinomial::dpbinom(i, p) <= d * rel.err)
        PoissonBinomial::ppbinom(x, p) +
          PoissonBinomial::ppbinom(length(p) - y, p, lower.tail = FALSE)
      } else {
        i <- seq.int(from = 0, to = floor(m))
        y <- sum(PoissonBinomial::dpbinom(i, p) <= d * rel.err)
        PoissonBinomial::ppbinom(y - 1, p) +
          PoissonBinomial::ppbinom(x - 1, p, lower.tail = FALSE)
      }
    }
  )
  p.value
}
