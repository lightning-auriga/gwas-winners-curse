#' @title
#' Compute asymptotic 95% confidence interval around debiased regression
#' coefficient estimate.
#'
#' @description
#' This just computes the 95% confidence interval. The original implementation
#' had a wildly overcomplicated zero solving method that is entirely
#' unnecessary given normal theory.
#'
#' @param beta.debiased Numeric parameter estimate.
#' @param stderr.debiased Numeric standard error of the mean.
#' @return Numeric vector, first entry lower 95% confidence bound,
#' second entry upper 95% confidence bound.
calculate.ci <- function(beta.debiased, stderr.debiased, p.thresh) {
  c(
    beta.debiased + qnorm(0.025) * stderr.debiased,
    beta.debiased + qnorm(0.975) * stderr.debiased
  )
}
