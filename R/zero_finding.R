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
#' @param p.thresh Numeric p-value boundary for signed Z statistic.
#' @return Numeric vector corresponding to provided p-value boundary.
calculate.ci <- function(beta.debiased, stderr.debiased, p.thresh) {
  beta.debiased + qnorm(p.thresh) * stderr.debiased
}

#' @title
#' Adjust mean of phenotype distribution to compensate for genetic effect.
#'
#' @description
#' The below equation is valid for linear regression only. The original
#' implementation tried to handle other GLM linkers, but that wasn't
#' truly supported in the rest of the package, so it is removed for
#' the time being.
#'
#' @param trait.characteristic Numeric mean of phenotype.
#' @param freq Numeric allele frequency of variant.
#' @param beta Numeric regression coefficient.
#' @return Numeric trait mean with expected genetic effect removed.
calculate.adjusted.trait <- function(trait.characteristic, freq, beta) {
  trait.characteristic - 2 * freq * beta
}

#' @title
#' Compute optimization function for iteratively debiasing regression coefficient.
#'
#' @description
#' This function was substantially more complicated in the original implementation,
#' but as with many other functions, this only actually supported linear regression
#' in practice, so other GLM linker support is removed for the time being.
#'
#' @param beta.debiased.init Numeric input attempted debiased value of regression coefficient.
#' @param p.thresh Numeric discovery p-value threshold.
#' @param beta.biased Numeric input biased regression coefficient.
#' @param stderr.biased Numeric input (possibly) biased standard error.
#' @return Numeric optimization function value for beta debiasing routine.
#' @importFrom stats dnorm pnorm qnorm
debiasing.func <- function(beta.debiased.init, p.thresh,
                           beta.biased, stderr.biased) {
  ## for linear regression, the standard error is invariant with respect to
  ## regression coefficient. so skip the entire concept of standard error
  ## rescaling until other GLM linker support is attempted.
  beta.debiased <- ifelse(abs(beta.debiased.init) < 1e-16, sign(beta.debiased.init) * 1e-16, beta.debiased.init)
  adjusted.stderr <- stderr.biased
  q <- beta.debiased / adjusted.stderr - qnorm(1.0 - p.thresh / 2.0)
  r <- -beta.debiased / adjusted.stderr - qnorm(1.0 - p.thresh / 2.0)
  beta.debiased + adjusted.stderr * (dnorm(q) - dnorm(r)) / (pnorm(q) + pnorm(r)) - beta.biased
}

#' @title
#' Compute MLE debiased regression coefficient.
#'
#' @description
#' Adjusts input regression coefficient according to Winner's Curse
#' logic and threshold used for discovery. Most of the logic present
#' in the original implementation has been removed to improve maintainability
#' and legibility.
#'
#' @details
#' Some of the original implementation was designed to support
#' non-linear GLMs, but the implementation was never finished and
#' left in an unstable state. That content has been removed in favor
#' of making the linear regression content transparent and testable.
#' The additional content may be added back in at a later date, maybe.
#'
#' Note that papers often misreport the actual p-value used to select
#' variants for replication. In some cases, a reasonable approximation
#' of the true threshold used can be derived from the observed p-values,
#' but if the investigators just spiked in manually selected variants,
#' the entire model somewhat breaks down.
#'
#' @param beta.biased Numeric input biased regression coefficient.
#' @param stderr.biased Numeric input coefficient standard error.
#' @param freq Numeric variant allele frequency.
#' @param trait.characteristic Numeric phenotype mean.
#' @param p.threshold.init Numeric reported p-value threshold used
#' to bring variants into replication.
#' @return Numeric maximum likelihood estimate of debiased regression coefficient.
#' @importFrom stats uniroot
debias.beta <- function(beta.biased,
                        stderr.biased,
                        freq,
                        trait.characteristic,
                        p.threshold.init) {
  beta.nosign <- abs(beta.biased)
  p.actual <- 2.0 * (1.0 - pnorm(beta.nosign / stderr.biased))
  p.threshold <- 0
  ## deal with the fact that paper-reported thresholds for inclusion in replication
  ## are often greatly exaggerated, and if included variants do not in fact pass the threshold,
  ## the calculation needs adjustment to not explode.
  if (p.actual > p.threshold.init) {
    p.threshold <- p.actual
    warning("pair ", beta.biased, " (", stderr.biased, ") corresponds to p = ",
      p.actual, " > ", p.threshold.init, ", adjusting accordingly",
      sep = ""
    )
  } else {
    p.threshold <- p.threshold.init
  }

  x.lo <- -beta.nosign / 4
  x.hi <- beta.nosign
  beta0 <- calculate.adjusted.trait(beta.nosign, freq, trait.characteristic)

  f.solver <- function(x) {
    debiasing.func(x, p.threshold, beta.nosign, stderr.biased)
  }
  beta.mle <- uniroot(f.solver, c(x.lo, x.hi))$root * sign(beta.biased)
  beta.mle
}
