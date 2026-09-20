## The live fitting path, on 400 synthetic rows (see helper-syntheticData.R). `ignitions` is a
## function of CMDsm and youngAge only, and `escapes` of lightning only, so what the fitted models
## must achieve is known without running them: they separate the outcomes (AUC = 1).

fitCache <- function() {
  withr::local_options(list(reproducible.cachePath = file.path(testPaths$cachePath, "runXGBOOST")),
                       .local_envir = parent.frame())
}

fitQuietly <- function(...) {
  res <- NULL
  utils::capture.output(suppressWarnings(suppressMessages(res <- runXGBOOST(...))))
  res
}

test_that("runXGBOOST returns one model per fold, fitted on the covariates only", {
  fitCache()
  dat <- scaleCovariates(makeIgnitionCovariates())
  set.seed(42)
  m <- fitQuietly(data.table::copy(dat), dig = "test", type = "ignition", nFolds = 5)

  expect_identical(names(m), c(paste0("Fold", 1:5), "rocs"))
  expect_identical(names(m$rocs), paste0("Fold", 1:5))
  for (fold in paste0("Fold", 1:5)) {
    ## pixelID, year and the response do not reach the model; the rest do, in alphabetical order
    vn <- stats::variable.names(m[[fold]])
    expect_setequal(vn, c("CMDsm", "lightning", "Pice_mar", "youngAge"))
    expect_identical(vn, sort(vn))
    pars <- attr(m[[fold]], "params")
    expect_identical(pars$objective, "reg:tweedie")
    expect_equal(pars[c("max_depth", "learning_rate", "reg_lambda")],
                 list(max_depth = 2, learning_rate = 0.15, reg_lambda = 0.5))
    expect_identical(nrow(attr(m[[fold]], "evaluation_log")), 100L)      # nrounds
  }
})

test_that("the ignition models recover the known signal", {
  fitCache()
  raw <- makeIgnitionCovariates()
  dat <- scaleCovariates(raw)
  set.seed(42)
  m <- fitQuietly(data.table::copy(dat), dig = "test", type = "ignition", nFolds = 5)

  ## the outcome is a function of two covariates, so every held-out fold is ranked perfectly
  expect_equal(unname(aucPerFold(m$rocs)), rep(1, 5))
  ## each ROC curve is on its fold's held-out rows: 400 rows in 5 folds, 61 of them positive
  nHeld <- vapply(m$rocs, function(r) length(r$response), numeric(1))
  expect_identical(sum(nHeld), 400)
  expect_true(all(nHeld %in% 79:81))
  expect_identical(sum(vapply(m$rocs, function(r) sum(r$response), numeric(1))), 61)

  ## predicted counts track the observed 0, 1 and 2 ignitions
  pred <- predict(m$Fold1, as.data.frame(dat)[, c("CMDsm", "lightning", "Pice_mar", "youngAge")])
  means <- tapply(pred, raw$ignitions, mean)
  expect_lt(means[["0"]], 0.01)
  expect_equal(means[["1"]], 1, tolerance = 0.15)
  expect_equal(means[["2"]], 2, tolerance = 0.05)
  ## value from origin/development at e3aa2a7 (xgboost 3.2.1.1)
  expect_equal(as.numeric(means), c(0.00042, 0.91199, 2.01936), tolerance = 1e-3)
})

test_that("runXGBOOST gives the same folds and models when repeated, and restores the RNG", {
  fitCache()
  dat <- scaleCovariates(makeIgnitionCovariates())
  set.seed(42)
  m1 <- fitQuietly(data.table::copy(dat), dig = "test", type = "ignition", nFolds = 5)
  afterFit <- stats::runif(1)
  set.seed(42)
  expect_identical(afterFit, stats::runif(1))          # the caller's random stream is untouched

  set.seed(999)                                        # a different seed outside: same folds inside
  m2 <- fitQuietly(data.table::copy(dat), dig = "test", type = "ignition", nFolds = 5)
  for (fold in paste0("Fold", 1:5))
    expect_identical(as.numeric(m1$rocs[[fold]]$predictor), as.numeric(m2$rocs[[fold]]$predictor))

  ## the folds are those of cvFolds() under the fixed seed 12345
  set.seed(12345)
  folds <- cvFolds(dat$ignitions, 5)
  for (fold in paste0("Fold", 1:5))
    expect_identical(as.numeric(m1$rocs[[fold]]$response),
                     as.numeric(pmin(1L, dat$ignitions[folds[[fold]]$keepEval])))
})

test_that("the number of folds", {
  fitCache()
  dat <- scaleCovariates(makeIgnitionCovariates())
  set.seed(42)
  m3 <- fitQuietly(data.table::copy(dat), dig = "test", type = "ignition", nFolds = 3)
  expect_identical(names(m3), c("Fold1", "Fold2", "Fold3", "rocs"))
  heldOut <- vapply(m3$rocs, function(r) length(r$response), numeric(1))
  expect_identical(sum(heldOut), 400)                  # 3 folds still cover every row once
})

test_that("the escape model uses only rows with an ignition, and `escapes` as response", {
  fitCache()
  raw <- makeIgnitionCovariates(escapes = TRUE)
  dat <- scaleCovariates(raw)
  set.seed(42)
  m <- fitQuietly(data.table::copy(dat), dig = "test", type = "escape", nFolds = 5)

  ## 61 rows have an ignition; 33 of them escaped (lightning >= 6)
  expect_identical(sum(raw$ignitions > 0), 61L)
  expect_identical(sum(raw$escapes), 33L)
  expect_identical(sum(vapply(m$rocs, function(r) length(r$response), numeric(1))), 61)
  expect_identical(sum(vapply(m$rocs, function(r) sum(r$response), numeric(1))), 33)
  ## neither response is a predictor
  expect_setequal(stats::variable.names(m$Fold1), c("CMDsm", "lightning", "Pice_mar", "youngAge"))
  ## escapes depend on lightning alone, so the held-out rows are ranked perfectly
  expect_equal(unname(aucPerFold(m$rocs)), rep(1, 5))

  withIgnition <- as.data.frame(dat)[raw$ignitions > 0, c("CMDsm", "lightning", "Pice_mar", "youngAge")]
  pred <- predict(m$Fold1, withIgnition)
  means <- tapply(pred, raw$escapes[raw$ignitions > 0], mean)
  expect_lt(means[["0"]], 0.05)
  expect_equal(means[["1"]], 1, tolerance = 0.1)
})

test_that("too few positives stops the fit before any model is trained", {
  fitCache()
  dat <- scaleCovariates(makeIgnitionCovariates())
  data.table::set(dat, NULL, "ignitions", 0L)
  data.table::set(dat, 5L, "ignitions", 1L)            # a single ignition
  set.seed(42)
  expect_error(fitQuietly(data.table::copy(dat), dig = "test", type = "ignition", nFolds = 5),
               "Too few ignitions to fit 5 cross-validation folds: 1 positive observation(s), and 1 fold(s)",
               fixed = TRUE)
})

test_that("buildModel with xgboost fits runXGBOOST on integer responses", {
  fitCache()
  dat <- scaleCovariates(makeIgnitionCovariates())
  data.table::set(dat, NULL, "ignitions", as.numeric(dat$ignitions))
  set.seed(42)
  res <- NULL
  utils::capture.output(suppressWarnings(suppressMessages(
    res <- buildModel(covariates = dat, type = "ignition",
                      dig = "test", modelAlgorithm = "xgboost", nFolds = 4)
  )))
  expect_identical(names(res), c(paste0("Fold", 1:4), "rocs"))            # nFolds is passed on
  expect_equal(unname(aucPerFold(res$rocs)), rep(1, 4))
  expect_type(dat$ignitions, "integer")                                  # converted by reference
  expect_identical(sum(dat$ignitions), 83L)                              # 39 * 1 + 22 * 2
})

test_that("buildModel refuses a non-xgboost algorithm (that path was removed)", {
  expect_error(buildModel(covariates = data.table::data.table(ignitions = 0L), type = "ignition",
                          dig = "test", modelAlgorithm = "glmmadaptive"),
               "non-xgboost path was removed", fixed = TRUE)
})
