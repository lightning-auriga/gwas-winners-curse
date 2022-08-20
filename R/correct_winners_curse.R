#' @title
#' Correct Winner's Curse in GWAS association data.
#'
#' @description
#' This function accepts a file containing summary information
#' describing the discovery association metrics for a set
#' of significant GWAS variants. The function will attempt
#' to adjust those values to correct for the apparent
#' "Winner's Curse" effect on those estimates, which depends
#' most notably on the p-value threshold used in discovery
#' to select variants for replication. Both MLE and MSE
#' corrections are emitted; generally, MLE has performed better
#' in testing, but that is not guaranteed to be accurate.
#'
#' @details
#' This function is a reimplementation of the original C++
#' program that accompanied the paper at
#' \url{https://doi.org/10.1371/journal.pgen.1006916}. The implementation
#' is largely improved and streamlined. As a brief summary of changes:
#'
#' - obviously, it's an R package now
#' - linear regression is explicitly the only method supported. that
#'   isn't actually a change from the original, as the cooked-in
#'   alternate GLM support was not directly part of the paper
#' - the input and output formats are slightly changed, which is fine
#'   because they weren't really document to begin with whoops.
#'   the formats are as follows:
#' - Input:
#'   - beta (discovery)
#'   - standard error (discovery)
#'   - sample size (discovery)
#'   - allele frequency (discovery)
#'   - phenotype distribution mean in study (optional)
#'   - p-value threshold to select variant for replication (optional, though
#'     required if phenotype distribution mean is specified as a column)
#'   - beta (replication)
#'   - standard error (replication)
#'   - sample size (replication)
#'   - allele frequency (replication)
#' - Output:
#'   - The same first columns as input. A header is added or slightly adjusted as needed.
#'     If no per-row trait mean/p-value threshold was provided in the input, they will be
#'     appended as per-row values in columns 9 and 10 of the output
#'   - beta (MLE adjustment for Winner's Curse)
#'   - 95% CI on beta, lower bound (MLE adjustment for Winner's Curse)
#'   - 95% CI on beta, upper bound (MLE adjustment for Winner's Curse)
#'   - beta (MSE adjustment for Winner's Curse)
#'   - 95% CI on beta, lower bound (MSE adjustment for Winner's Curse)
#'   - 95% CI on beta, upper bound (MSE adjustment for Winner's Curse)
#'
#' - Additional modifications to input and output format are available based
#'   on formal parameters to the function, as described in the parameter documentation.
#' - In the original package, standard error of the input regression coefficient was being
#'   reestimated. This was originally intended in advance of support for non-linear GLM,
#'   but that support never manifested, and as such, the correction was actually somewhat
#'   counterproductive, seeing as:
#'   - there is a simple, fast closed-form solution to standard error
#'     for additive genetic effect under HWE for linear models (so it just went slow for no reason)
#'   - any deviation from HWE (which is expected for some significant variants) caused
#'     the correction to be inaccurate
#'   - linear model genetic effect standard error is independent of beta (yes, it's true), and
#'     as such, no correction is actually needed (again, only for linear regression)
#'   In light of the above, the correction method is wholly removed, and will only be added
#'   back in if other GLM support is attempted at some later date.
#'
#' @param input.file Character vector name of input data file.
#' @param output.file Character vector name of output file
#' to which results should be reported.
#' @param trait.mean Numeric mean of trait used for discovery.
#' If NA, the function will expect mean values to be specified in
#' each row of the input data file.
#' @param discovery.threshold Numeric threshold on p-value used
#' in discovery phase to select variants for replication. If NA,
#' the function will expect threshold values to be specified in
#' each row of the input data file.
#' @param header Logical indicating whether input data file
#' has a header.
#' @param sep Character vector indicating delimiter in input data.
#' @return NULL
#' @importFrom utils read.table write.table
#' @export
correct.winners.curse <- function(input.file,
                                  output.file,
                                  trait.mean = 0,
                                  discovery.threshold = 5e-8,
                                  header = TRUE,
                                  sep = "\t") {
  stopifnot(length(input.file) == 1)
  stopifnot(is.character(input.file))
  stopifnot(file.exists(input.file))
  stopifnot(length(output.file) == 1)
  stopifnot(is.character(output.file))
  stopifnot(length(trait.mean) == 1)
  stopifnot(is.numeric(trait.mean) || is.na(trait.mean))
  stopifnot(length(discovery.threshold) == 1)
  stopifnot(is.numeric(discovery.threshold) || is.na(discovery.threshold))

  ## load input data
  data <- read.table(input.file, header = header, sep = sep, stringsAsFactors = FALSE)
  stopifnot(nrow(data) > 0)
  stopifnot(ncol(data) %in% c(8, 10))
  if (ncol(data) == 10) {
    colnames(data) <- c(
      "discovery.beta", "discovery.se", "discovery.n", "discovery.freq",
      "trait.mean", "discovery.threshold",
      "replication.beta", "replication.se", "replication.n", "replication.freq"
    )
  } else {
    colnames(data) <- c(
      "discovery.beta", "discovery.se", "discovery.n", "discovery.freq",
      "replication.beta", "replication.se", "replication.n", "replication.freq"
    )
    data$trait.mean <- trait.mean
    data$discovery.threshold <- discovery.threshold
  }
  ## debias beta
  data$debiased.beta.mle <- sapply(seq_len(nrow(data)), function(i) {
    debias.beta(
      data[i, "discovery.beta"],
      data[i, "discovery.se"],
      data[i, "discovery.freq"],
      data[i, "trait.mean"],
      data[i, "discovery.threshold"]
    )
  })
  ## compute 95% confidence interval around debiased beta
  data$l95.mle <- calculate.ci(data[, "debiased.beta.mle"], data[, "discovery.se"], 0.025)
  data$u95.mle <- calculate.ci(data[, "debiased.beta.mle"], data[, "discovery.se"], 0.975)
  ## compute MSE variants of debiased results
  data$debiased.beta.mse <- compute.beta.mse(
    data[, "discovery.beta"],
    data[, "discovery.se"],
    data[, "debiased.beta.mle"]
  )
  data$l95.mse <- compute.ci.mse(
    data[, "discovery.beta"],
    data[, "discovery.se"],
    data[, "l95.mle"],
    0.025
  )
  data$u95.mse <- compute.ci.mse(
    data[, "discovery.beta"],
    data[, "discovery.se"],
    data[, "u95.mle"],
    0.975
  )
  ## fix precision in computed variables
  for (colname in c(
    "debiased.beta.mle", "l95.mle", "u95.mle",
    "debiased.beta.mse", "l95.mse", "u95.mse"
  )) {
    data[, colname] <- signif(data[, colname], 5)
  }
  ## emit results
  write.table(data, output.file,
    row.names = FALSE, col.names = TRUE, quote = FALSE,
    sep = ifelse(grepl("\\.csv$", output.file), ",", "\t")
  )
}
