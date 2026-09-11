## ELF 3.1.2 had a single ignition in 2002-2022. caret::createFolds() put it in one validation fold;
## xgboost trains each fold's model on the rows outside its eval_set, so that fold's model saw only
## zeros and predicted NaN for every row. pROC::roc() then dropped every row and stopped with
## "No control observation", which the single-outcome guard in rocPerFold() does not catch.

sourceModule <- function() {
  e <- new.env(parent = globalenv())
  e$defineModule <- function(...) invisible(NULL)
  e$defineParameter <- e$expectsInput <- e$createsOutput <- function(...) NULL
  suppressWarnings(sys.source(
    testthat::test_path("..", "..", "fireSense_IgnitionFit.R"), envir = e, keep.source = FALSE))
  e
}

test_that("xgboost's tweedie model predicts NaN when it trains on zeros only", {
  skip_if_not_installed("xgboost")
  x <- data.frame(a = seq(0, 1, length.out = 50))
  m <- xgboost::xgboost(x = x, y = rep(0, 50), objective = "reg:tweedie", nrounds = 5,
                        verbosity = 0, nthread = 1)
  expect_true(all(is.nan(predict(m, x))))
})

test_that("cvFolds spreads the positives so every fold trains with some", {
  skip_if_not_installed("caret")
  e <- sourceModule()
  for (positives in c(2L, 3L, 7L)) {
    response <- c(rep(3, positives), rep(0, 2000)) # counts, not only 0/1
    for (seed in 1:20) {
      set.seed(seed)
      folds <- e$cvFolds(response, nFolds = 5)
      expect_length(folds, 5)
      trainPositives <- vapply(folds, function(f) sum(response[setdiff(f$keepAll, f$keepEval)] > 0), numeric(1))
      expect_true(all(trainPositives >= 1))
      heldOutPositives <- vapply(folds, function(f) sum(response[f$keepEval] > 0), numeric(1))
      expect_lte(diff(range(heldOutPositives)), 1)
      expect_no_error(e$stopIfFoldsLackPositives(response, folds, "ignition"))
    }
  }
})

test_that("a fold that would train without any positive stops the fit with a clear message", {
  e <- sourceModule()
  fold <- function(keepEval) list(keepAll = 1:10, keepEval = keepEval)
  folds <- list(fold(1:2), fold(3:10))

  ## one positive, in the first fold's held-out rows: that fold trains on zeros only
  expect_error(e$stopIfFoldsLackPositives(c(1, rep(0, 9)), folds, "ignition"),
               "Too few ignitions to fit 2 cross-validation folds: 1 positive observation")

  ## a positive outside each fold's held-out rows: nothing to stop
  expect_no_error(e$stopIfFoldsLackPositives(c(1, rep(0, 8), 1), folds, "ignition"))

  ## the training rows are keepAll without keepEval, as for time-ordered folds
  expect_error(e$stopIfFoldsLackPositives(c(0, 0, 1, rep(0, 7)),
                                          list(list(keepAll = 1:4, keepEval = 3:4)), "escape"),
               "Too few escapes")
})
