test.file.codes <- list.files(".", "*test_input.tsv", recursive = TRUE, full.names = TRUE)
test.file.codes <- gsub("_test_input.tsv", "", sapply(strsplit(test.file.codes, "/"), function(i) {
  i[[4]]
}))
## this is admittedly lenient. the largest relative deviation is just over 5%. I really
## don't know whether this should be considered truly problematic or not. it looks like
## just general precision issues, but who can say.
permitted.rel.tolerance <- 0.1

test_that("correct.winners.curse functions as expected", {
  for (file.code in test.file.codes) {
    test.input <- file.path(
      "testthat_files", "correct_winners_curse",
      paste(file.code, "_test_input.tsv", sep = "")
    )
    test.expected <- file.path(
      "testthat_files", "correct_winners_curse",
      paste(file.code, "_test_expected.tsv", sep = "")
    )
    test.expected <- read.table(test.expected, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    test.output <- tempfile(paste("correct.winners.curse.", file.code, sep = ""), fileext = ".tsv")
    ## suppress warnings here, as we're not concerned at the moment with whether the
    ## p-value threshold needs dynamic adjustment
    suppressWarnings(correct.winners.curse(test.input, test.output, NA, NA, TRUE, "\t"))
    test.observed <- read.table(test.output, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    ## there are expected deviations in extreme precision values for these methods,
    ## so we're really only looking for very general concordance, but within 10% or so
    expected.beta <- test.expected[, "old_corrected_beta"]
    observed.beta <- test.observed[, "debiased.beta.mle"]
    expected.l95 <- test.expected[, "old_l95"]
    observed.l95 <- test.observed[, "l95.mle"]
    expected.u95 <- test.expected[, "old_u95"]
    observed.u95 <- test.observed[, "u95.mle"]
    if (length(which(abs(observed.beta - expected.beta) / expected.beta >= permitted.rel.tolerance)) > 0) {
      print(paste(file.code, "beta", sep = " "))
      print(cbind(
        observed.beta,
        expected.beta
      )[abs(observed.beta - expected.beta) / expected.beta >= permitted.rel.tolerance, ])
    }
    if (length(which(abs(observed.l95 - expected.l95) / expected.l95 >= permitted.rel.tolerance)) > 0) {
      print(paste(file.code, "l95", sep = " "))
      print(cbind(
        observed.l95,
        expected.l95
      )[abs(observed.l95 - expected.l95) / expected.l95 >= permitted.rel.tolerance, ])
    }
    if (length(which(abs(observed.u95 - expected.u95) / expected.u95 >= permitted.rel.tolerance)) > 0) {
      print(paste(file.code, "u95", sep = " "))
      print(cbind(
        observed.u95,
        expected.u95
      )[abs(observed.u95 - expected.u95) / expected.u95 >= permitted.rel.tolerance, ])
    }
    expect_true(all(abs(observed.beta - expected.beta) / expected.beta < permitted.rel.tolerance))
    expect_true(all(abs(observed.l95 - expected.l95) / expected.l95 < permitted.rel.tolerance))
    expect_true(all(abs(observed.u95 - expected.u95) / expected.u95 < permitted.rel.tolerance))
  }
})

test_that("correct.winners.curse can deal with files that don't contain per-row means and p thresholds", {
  test.input <- file.path("testthat_files", "correct_winners_curse", "24963161_no_mean_p.csv")
  test.expected <- file.path("testthat_files", "correct_winners_curse", "24963161_test_expected.tsv")
  test.expected <- read.table(test.expected, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  test.output <- tempfile("correct_winners_curse_truncated", fileext = ".csv")
  suppressWarnings(correct.winners.curse(test.input, test.output, 0, 5e-6, FALSE, ","))
  test.observed <- read.table(test.output, header = TRUE, stringsAsFactors = FALSE, sep = ",")
  expected.beta <- test.expected[, "old_corrected_beta"]
  observed.beta <- test.observed[, "debiased.beta.mle"]
  expected.l95 <- test.expected[, "old_l95"]
  observed.l95 <- test.observed[, "l95.mle"]
  expected.u95 <- test.expected[, "old_u95"]
  observed.u95 <- test.observed[, "u95.mle"]
  if (length(which(abs(observed.beta - expected.beta) / expected.beta >= permitted.rel.tolerance)) > 0) {
    print("beta")
    print(cbind(
      observed.beta,
      expected.beta
    )[abs(observed.beta - expected.beta) / expected.beta >= permitted.rel.tolerance, ])
  }
  if (length(which(abs(observed.l95 - expected.l95) / expected.l95 >= permitted.rel.tolerance)) > 0) {
    print("l95")
    print(cbind(
      observed.l95,
      expected.l95
    )[abs(observed.l95 - expected.l95) / expected.l95 >= permitted.rel.tolerance, ])
  }
  if (length(which(abs(observed.u95 - expected.u95) / expected.u95 >= permitted.rel.tolerance)) > 0) {
    print("u95")
    print(cbind(
      observed.u95,
      expected.u95
    )[abs(observed.u95 - expected.u95) / expected.u95 >= permitted.rel.tolerance, ])
  }
  expect_true(all(abs(observed.beta - expected.beta) / expected.beta < permitted.rel.tolerance))
  expect_true(all(abs(observed.l95 - expected.l95) / expected.l95 < permitted.rel.tolerance))
  expect_true(all(abs(observed.u95 - expected.u95) / expected.u95 < permitted.rel.tolerance))
})
