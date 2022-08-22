#' @title
#' Compute approximate phenotypic variance explained by a genetic variant.
#'
#' @description
#' Given certain very important assumptions, notably that the genetic effect
#' of a variant is nearly 0, the genetic variance explained under additivity
#' for biallelic variants is simply 2B^2f(1-f), for B regression coefficient,
#' f allele frequency.
#'
#' @details
#' Note that the genetic variance explained described above breaks down when
#' the regression coefficient substantially exceeds 0. It is slightly ridiculous
#' to consider an equation that self-destructs when it actually has something
#' to estimate, but there you go.
#'
#' @param beta Numeric additive regression coefficient for a biallelic variant.
#' @param freq Numeric allele frequency for a biallelic variant.
#' @return Numeric approximate variance explained for the variant, assuming the
#' association is small in magnitude.
#' @export
#' @seealso compute.power
#' @examples
#' freq <- 0.25
#' beta <- 0.05
#' varexp <- compute.variance.explained(beta, freq)
compute.variance.explained <- function(beta, freq) {
  2 * beta^2 * freq * (1 - freq)
}

#' @title
#' Compute power to detect a signal given certain assumptions.
#'
#' @description
#' The power to detect an association is, in short, the probability
#' of the variant exceeding the required test statistic threshold
#' as set under the null hypothesis when the alternative hypothesis
#' is true. In formal terms, this becomes 1 - \chi^2_1((\chi^2_1)^{-1}(1 - \alpha), ncp_v)
#' where the noncentrality parameter is N\frac{variance explained}{residual variance}.
#' How fine grained you want to get with this calculation is really up to personal
#' preference.
#'
#' @details
#' Which version of these values you provide to the function depends on the question
#' you want to ask. In the context of predicting replication, you probably want
#' to provide the _debiased discovery beta_ and the _replication sample size and frequency_.
#' But different permutations of these values will answer different questions.
#'
#' @param beta Numeric additive regression coefficient for a biallelic variant. In the
#' context of power calculations, this is meant to be a proxy for the true variant effect.
#' @param freq Numeric allele frequency for a biallelic variant. This should probably
#' be the frequency in the population (or sample) you'll be using to try to detect a signal,
#' and should be as accurate and precise as possible.
#' @param sample.size Numeric number of subjects in study sample. This should be as accurate
#' as possible; so for example, in the context of GWAS meta-analysis, do not simply provide
#' the _maximum_ sample size available in theory, but rather the actual sample size available
#' for the variant in question.
#' @param p.threshold Numeric p-value of significance required to call a variant
#' as detected. Depending on what phase of a study you're in, this value may differ substantially.
#' For studies testing millions of variants in discovery, this will likely be 5*10^-8 (genome-wide
#' significance) or lower; for replication, this may be more like \frac{0.05}{variants brought into replication}.
#' @return Numeric approximate power to detect a signal given the input descriptors of
#' that signal and study sample.
#' @export
#' @seealso compute.variance.explained
#' @examples
#' freq <- 0.25
#' beta <- 0.1
#' sample.size <- 10000
#' p.threshold <- 5e-8
#' compute.power(beta, freq, sample.size, p.threshold)
compute.power <- function(beta, freq, sample.size, p.threshold) {
  test.stat <- qchisq(1 - p.threshold, 1, lower.tail = TRUE)
  varexp <- compute.variance.explained(beta, freq)
  ncp <- sample.size * varexp / (1 - varexp)
  power.to.detect <- pchisq(test.stat, 1, lower.tail = FALSE, ncp = ncp)
  power.to.detect
}

#' @title
#' For a set of variants with data about both discovery and replication,
#' determine the number of variants expected to successfully replicate.
#'
#' @description
#' This function computes power to detect each variant in replication, and
#' combines this information across all variants to create an expected value
#' of number of variants to replicate.
#'
#' @details
#' This step is really challenging to interpret correctly. If all values are
#' correctly specified, the quality (or lack thereof) of this calculation will
#' be determined by the extent to which (1) the discovery regression coefficient
#' is biased by, for example, the Winner's Curse, and (2) the extent to which
#' shared genetic effect is actually present across the discovery and replication
#' samples.
#'
#' However, it is in fact most likely that the values will _not_ be correctly
#' specified in this step, in which case the interpretation is much more vague.
#' Specifically, the sample size and allele frequency data available for most studies
#' is both inaccurate and imprecise, due to a variety of challenging factors.
#' Any amount of misspecification of these values will cause problems; but very notably,
#' the tendency of meta-analyses (or the GWAS catalog) to specify _maximum_ sample
#' size instead of _actual per-variant_ sample size causes systematic overestimation
#' of the power to replicate. There is no easy answer to this problem, and unfortunately
#' resources like the GWAS catalog are actually moving farther away from any rigorous
#' standard of data collection in favor of a simplified standardization of input information,
#' to maximize collection. Thus it seems likely that it will actually become
#' less possible over time to get systematic replication assessments in the field,
#' and interested investigators should conduct these calculations within their own
#' datasets where they have access to improved information.
#'
#' @param variant.effect Numeric variant effect size. In theory, this is meant
#' to be the true underlying effect estimate, though in practice this will likely
#' be estimated from another dataset.
#' @param study.frequency Numeric allele frequency of variant within the study sample.
#' Deviations of this number from the actual study allele frequency will negatively
#' impact the replication count estimation. In particular, reference allele frequency
#' from external datasets is discouraged for use here.
#' @param study.sample.size Numeric number of subjects in study sample.
#' Deviations of this number from the actual study sample size will negatively
#' impact the replication count estimation. In particular, this should _not_
#' be _maximum sample size_ but rather _actual sample size_ per variant.
#' @return Numeric expected number of variants from the combined variant dataset
#' that are expected to replicate, given the provided study parameters.
#' @export
#' @seealso compute.power, compute.variance.explained
#' @examples
#' beta <- c(0.1, 0.2)
#' f <- c(0.22, 0.12)
#' n <- c(10000, 5000)
#' p.thresh <- c(1e-5, 1e-6)
#' expected.replication.count(beta, f, n, p.thresh)
expected.replication.count <- function(variant.effect,
                                       study.frequency,
                                       study.sample.size,
                                       p.threshold) {
  ## expected replication count is merely the sum of the power.
  ## a more complicated question is whether this outcome was likely
  ## given the sequence of power estimates.
  power.to.detect <- compute.power(
    variant.effect,
    study.frequency,
    study.sample.size,
    p.threshold
  )
  expected.count <- sum(power.to.detect)
  expected.count
}
