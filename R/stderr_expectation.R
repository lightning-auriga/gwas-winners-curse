#' @title
#' Compute expected standard error for genetic additive
#' coefficient in linear regression.
#'
#' @description
#' The linear regression standard error for the additive
#' model is independent of the regression coefficient,
#' and has a simple closed-form solution.
#'
#' @details
#' The original implementation of this function had theoretical
#' support for other GLM linkers, but the rest of the package
#' did not support that, so I'm not really sure what it was
#' doing here to begin with.
#'
#' This equation assumes Hardy Weinberg, which basically
#' just renders it a worse approximation of the error than
#' the input error provided by the user. Nevertheless, this
#' is recorded here for posterity.
#'
#' @param freq Numeric allele frequency.
#' @param n.init Integer input sample size.
#' @return Numeric predicted standard error for
#' linear regression additive genetic effect estimate.
linear.stderr.expectation <- function(freq, n.init) {
  1 / sqrt(n.init * 2 * freq * (1 - freq))
}
