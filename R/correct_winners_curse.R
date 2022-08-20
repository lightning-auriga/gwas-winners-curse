#' @title
#' Correct Winner's Curse in GWAS association data
#'
#' @description
#' TBD
#'
#' @param input.file Character vector name of input data file.
#' @param output.file Character vector name of output file
#' to which results should be reported.
#' @param regression.type Character vector name of regression
#' type used to generate association data.
#' @param trait.mean Numeric mean of trait used for discovery.
#' If NA, the function will expect mean values to be specified in
#' each row of the input data file.
#' @param stderr.precision Numeric requested precision of standard
#' error calculation.
#' @param discovery.threshold Numeric threshold on p-value used
#' in discovery phase to select variants for replication.
#' @param header Logical indicating whether input data file
#' has a header.
#' @param invariant.stderr Logical indicating whether to use
#' the simplifying assumption that the standard error does not
#' change along with shrinking effect estimate.
#' @return NULL
#' @export
correct.winners.curse <- function(input.file,
                                  output.file,
                                  regression.type = c("linear", "logistic", "poisson"),
                                  trait.mean = NA,
                                  stderr.precision = NA,
                                  discovery.threshold = 5e-8,
                                  header = TRUE,
                                  invariant.stderr = FALSE) {
  stopifnot(length(input.file) == 1)
  stopifnot(is.character(input.file))
  stopifnot(length(output.file) == 1)
  stopifnot(is.character(output.file))
  regression.type <- match.arg(regression.type)
}
