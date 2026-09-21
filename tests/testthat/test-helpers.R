## Small pure functions: values are worked out by hand in the comments.

test_that("functionNameHelper joins its arguments with '_' by default", {
  expect_identical(functionNameHelper("xgboost", "ignition", 3L), "xgboost_ignition_3")
  expect_identical(functionNameHelper("a", "b", sep = "-"), "a-b")
  ## an empty piece is kept, as in the escape plot filename `FuelClimate__predicted`
  expect_identical(functionNameHelper("FuelClimate", "", "predicted"), "FuelClimate__predicted")
  expect_identical(functionNameHelper("f", 1:2), c("f_1", "f_2"))
})

test_that("rocPerFold and aucPerFold give the AUC of each fold", {
  fold <- function(outcome, pred) list(valData = data.frame(ignitions = outcome, predTweedie = pred))
  ## cases predicted 0.2 and 0.4, controls 0.1 and 0.3: of the 4 case-control pairs, 3 have the
  ## case above the control (0.2 > 0.1, 0.4 > 0.1, 0.4 > 0.3), so AUC = 3/4.
  threeQuarters <- fold(c(0, 1, 0, 1), c(0.1, 0.2, 0.3, 0.4))
  ## counts above 1 are presences: the same fold with counts 3 and 2 has the same AUC
  counts <- fold(c(0, 3, 0, 2), c(0.1, 0.2, 0.3, 0.4))
  ## every case above every control: AUC = 1
  perfect <- fold(c(0, 0, 1, 1), c(0.1, 0.2, 0.8, 0.9))
  rocs <- suppressMessages(rocPerFold(list(a = threeQuarters, b = counts, c = perfect), "ignitions"))
  expect_identical(names(rocs), c("a", "b", "c"))
  expect_equal(aucPerFold(rocs), c(a = 0.75, b = 0.75, c = 1))
  expect_identical(as.integer(rocs$b$response), c(0L, 1L, 0L, 1L))
  expect_equal(as.numeric(rocs$a$predictor), c(0.1, 0.2, 0.3, 0.4))

  ## the response column is looked up by name
  esc <- list(valData = data.frame(ignitions = c(1, 1, 1, 1), escapes = c(0, 1, 0, 1),
                                   predTweedie = c(0.1, 0.2, 0.3, 0.4)))
  expect_equal(aucPerFold(suppressMessages(rocPerFold(list(esc), "escapes"))), 0.75)
})

test_that("a fold with one outcome is skipped with a message naming the value", {
  one <- list(valData = data.frame(escapes = c(2, 1, 1), predTweedie = c(0.2, 0.5, 0.9)))
  expect_message(rocs <- rocPerFold(list(one), "escapes"),
                 "a validation fold of escapes holds only the value 1; AUC is undefined", fixed = TRUE)
  expect_identical(rocs, list(NULL))
  expect_identical(aucPerFold(rocs), NA_real_)
  expect_identical(aucPerFold(list()), numeric(0))
})

test_that("meanRocMessage reports the mean over the folds that have an AUC", {
  expect_identical(meanRocMessage(c(0.75, 1)), "mean roc:  0.875")           # (0.75 + 1) / 2
  expect_identical(meanRocMessage(c(0.7, 0.8, 0.9)), "mean roc:  0.8")
  expect_identical(meanRocMessage(c(0.75, NA)), "mean roc:  0.75 (over 1 of 2 folds)")
  expect_identical(meanRocMessage(c(NA, 0.5, 1, NA)), "mean roc:  0.75 (over 2 of 4 folds)")
  expect_identical(meanRocMessage(c(NA_real_, NA_real_)),
                   "mean roc:  not computed (no validation fold had both outcomes)")
  expect_identical(meanRocMessage(2 / 3), "mean roc:  0.667")               # 3 significant digits
})

test_that("cvFolds partitions the rows and every fold keeps all rows in keepAll", {
  response <- c(rep(0, 40), rep(1, 6), rep(4, 4))        # 50 rows, 10 positive
  set.seed(7)
  folds <- cvFolds(response, nFolds = 5)
  expect_identical(names(folds), paste0("Fold", 1:5))
  for (f in folds) {
    expect_identical(names(f), c("keepAll", "keepEval"))
    expect_identical(f$keepAll, 1:50)
  }
  held <- unlist(lapply(folds, `[[`, "keepEval"), use.names = FALSE)
  expect_identical(sort(held), 1:50)                     # each row held out exactly once
  ## stratified on positive / not: 10 positives and 40 zeros over 5 folds is 2 and 8 in each
  expect_identical(unname(vapply(folds, function(f) sum(response[f$keepEval] > 0), numeric(1))),
                   rep(2, 5))
  expect_identical(unname(lengths(lapply(folds, `[[`, "keepEval"))), rep(10L, 5))

  ## the number of folds is honoured
  set.seed(7)
  expect_length(cvFolds(response, nFolds = 2), 2L)
})

test_that("stopIfFoldsLackPositives counts positives and failing folds in its message", {
  fold <- function(keepEval) list(keepAll = 1:6, keepEval = keepEval)
  response <- c(2, 1, 0, 0, 0, 0)
  ## fold 1 holds out both positives and trains on rows 3:6, all zero; folds 2 and 3 are fine
  expect_error(
    stopIfFoldsLackPositives(response, list(fold(1:2), fold(3:4), fold(5:6)), "escape"),
    "Too few escapes to fit 3 cross-validation folds: 2 positive observation(s), and 1 fold(s) would train without any.",
    fixed = TRUE)
  ## no positive at all: every fold fails
  expect_error(stopIfFoldsLackPositives(rep(0, 6), list(fold(1:3), fold(4:6)), "ignition"),
               "0 positive observation(s), and 2 fold(s)", fixed = TRUE)
  ## one positive left in the training rows of each fold is enough
  expect_null(stopIfFoldsLackPositives(response, list(fold(1L), fold(2L), fold(3:6)), "ignition"))
})
