#' @title
#' Correct Winner's Curse in GWAS association data.
#'
#' @description
#' TBD
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
