## A validation fold holding a single outcome must not abort the fit.
##
## `pROC::roc()` stops with "'response' must have two levels" when every observation has
## the same outcome. For escape in a small study area that is ordinary -- every fire
## escaped, or none did -- and it used to take down the entire ignition/escape fit from
## inside `purrr::pmap()`. On the run this was found on it cost a worker per attempt: the
## job failed, the worker parked at a prompt, the row went back to PENDING, and the next
## worker picked it up and parked too.
##
## The AUCs are a diagnostic. They are printed and returned beside the models; the fitted
## model (`m$mod`) does not depend on them. So a fold with one outcome reports NA.

sourceModule <- function() {
  e <- new.env(parent = globalenv())
  e$defineModule <- function(...) invisible(NULL)
  e$defineParameter <- e$expectsInput <- e$createsOutput <- function(...) NULL
  suppressWarnings(sys.source(
    testthat::test_path("..", "..", "fireSense_IgnitionFit.R"), envir = e, keep.source = FALSE))
  e
}

fold <- function(outcome, pred) list(valData = data.frame(escape = outcome, predTweedie = pred))

test_that("a single-outcome fold yields NA instead of stopping the fit", {
  skip_if_not_installed("pROC")
  e <- sourceModule()

  set.seed(1)
  mixed <- fold(c(0, 0, 1, 1, 0, 1), c(0.1, 0.2, 0.8, 0.9, 0.15, 0.7))
  allOnes  <- fold(rep(1, 6), runif(6))
  allZeros <- fold(rep(0, 6), runif(6))

  ## the failing shape: pROC alone stops here
  expect_error(pROC::roc(rep(1, 6), runif(6)), "two levels")

  ## ...and the module no longer does
  rocs <- suppressMessages(e$rocPerFold(list(mixed, allOnes, allZeros), "escape"))
  expect_length(rocs, 3L)
  expect_s3_class(rocs[[1]], "roc")
  expect_null(rocs[[2]])
  expect_null(rocs[[3]])

  aucs <- e$aucPerFold(rocs)
  expect_false(is.na(aucs[1]))
  expect_true(all(is.na(aucs[2:3])))

  ## the mean is over the folds that have one, and says so
  expect_match(e$meanRocMessage(aucs), "over 1 of 3 folds", fixed = TRUE)
  ## every fold single-outcome: no AUC at all, still no error
  expect_match(e$meanRocMessage(c(NA_real_, NA_real_)), "not computed", fixed = TRUE)
  ## the ordinary case is unchanged -- no parenthetical
  expect_false(grepl("of", e$meanRocMessage(c(0.8, 0.9)), fixed = TRUE))
})
