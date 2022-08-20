#' @title
#' Compute MSE adjustment for beta correction.
#'
#' @description
#' Computes weighted average of biased and strictly debiased
#' regression coefficient to report as intermediate result.
#'
#' @details
#' Note that this was apparently completely disabled in the original
#' implementation, such that the output reflected the MLE correction only.
#'
#' @param biased.beta Numeric input regression coefficient from discovery.
#' @param biased.se Numeric input coefficient standard error from discovery.
#' @param debiased.beta Numeric MLE adjusted regression coefficient.
#' @return Numeric weighted average coefficient reflecting MSE adjustment.
compute.beta.mse <- function(biased.beta,
                             biased.se,
                             debiased.beta) {
  k <- biased.se^2 / (biased.se^2 + (biased.beta - debiased.beta)^2)
  k * biased.beta + (1 - k) * debiased.beta
}

#' @title
#' Compute MSE adjustment for confidence interval correction.
#'
#' @description
#' Computes weighted average of biased and strictly debiased
#' regression coefficient confidence interval to report as intermediate result.
#'
#' @details
#' Note that this was apparently completely disabled in the original
#' implementation, such that the output reflected the MLE correction only.
#'
#' @param biased.beta Numeric input regression coefficient from discovery.
#' @param biased.se Numeric input coefficient standard error from discovery.
#' @param existing.ci Numeric MLE adjusted confidence interval.
#' @param p.thresh Numeric p-value boundary for signed Z statistic.
#' @return Numeric weighted average confidence interval reflecting MSE adjustment.
compute.ci.mse <- function(biased.beta,
                           biased.se,
                           existing.ci,
                           p.thresh) {
  biased.ci <- calculate.ci(biased.beta, biased.se, p.thresh)
  k <- biased.se^2 / (biased.se^2 + (biased.ci - existing.ci)^2)
  k * biased.ci + (1 - k) * existing.ci
}
