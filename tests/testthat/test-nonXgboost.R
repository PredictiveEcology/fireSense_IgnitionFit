## The non-xgboost code in buildModel() cannot be reached from the module today (see the PR), but
## it is still there, so its formula handling and model choice are pinned here. The model fitter is
## replaced by a Poisson glm that records the formulas it is given.

glmData <- function() {
  i <- 1:200
  data.table::data.table(
    pixelID = i, year = rep(2001:2004, 50), yearChar = as.character(rep(2001:2004, 50)),
    CMDsm = ((i * 37) %% 101) / 50,                     # 0 .. 2
    youngAge = as.numeric(i %% 2),
    Pice_mar = ((i * 53) %% 97) / 50,
    ignitions = 0
  )
}

recordingGlm <- function(record) {
  function(nam, form, dat, family, type, climVar) {
    record$forms[[nam]] <- form
    eval(bquote(stats::glm(.(form), data = dat, family = stats::poisson())))
  }
}

test_that("buildModel derives the no-interaction and climate-only formulas from the full one", {
  record <- new.env()
  testthat::local_mocked_bindings(runGLMMAdaptiveWithSimplifications = recordingGlm(record))
  dat <- glmData()
  data.table::set(dat, NULL, "ignitions", round(exp(dat$CMDsm * dat$Pice_mar)))

  best <- suppressMessages(
    buildModel(covariates = dat, formula = "ignitions ~ youngAge + CMDsm:Pice_mar + youngAge:CMDsm - 1",
               type = "ignition", climVar = "CMDsm", family = quote(poisson(link = "log")),
               digestOfData = list(), modelAlgorithm = "glmmadaptive"))

  expect_identical(names(record$forms), c("full", "InterceptOnly", "climateOnly"))
  expect_identical(messageFormulaFn(record$forms$full),
                   "ignitions ~ youngAge + CMDsm:Pice_mar + youngAge:CMDsm - 1")
  ## terms with ':' are dropped ...
  expect_identical(messageFormulaFn(record$forms$InterceptOnly), "ignitions ~ youngAge - 1")
  ## ... and the climate variable is added back as a main effect
  expect_identical(messageFormulaFn(record$forms$climateOnly), "ignitions ~ youngAge - 1 + CMDsm")

  ## the response is exp(CMDsm * Pice_mar), which only the full formula can express
  expect_identical(messageFormulaFn(stats::formula(best)),
                   "ignitions ~ youngAge + CMDsm:Pice_mar + youngAge:CMDsm - 1")
  expect_equal(unname(stats::coef(best)[["CMDsm:Pice_mar"]]), 1, tolerance = 0.02)
  ## the fitted model is returned without its data
  expect_null(best$y)
  expect_identical(ls(best$data), character(0))
  ## responses and year are made integer, yearChar a factor, by reference
  expect_type(dat$ignitions, "integer")
  expect_type(dat$year, "integer")
  expect_identical(levels(dat$yearChar), as.character(2001:2004))
})

test_that("buildModel keeps the simpler model when the interactions add nothing", {
  record <- new.env()
  testthat::local_mocked_bindings(runGLMMAdaptiveWithSimplifications = recordingGlm(record))
  dat <- glmData()
  data.table::set(dat, NULL, "ignitions", round(exp(1 + dat$CMDsm)))     # climate only

  best <- suppressMessages(
    buildModel(covariates = dat, formula = "ignitions ~ 1 + CMDsm:Pice_mar",
               type = "ignition", climVar = "CMDsm", family = quote(poisson(link = "log")),
               digestOfData = list(), modelAlgorithm = "glmmadaptive"))
  expect_identical(messageFormulaFn(stats::formula(best)), "ignitions ~ 1 + CMDsm")
  expect_equal(unname(stats::coef(best)), c(1, 1), tolerance = 0.02)     # log(mean) = 1 + 1 * CMDsm

})

test_that("an algorithm buildModel does not know stops with a message", {
  dat <- glmData()
  expect_error(
    buildModel(covariates = dat, formula = "ignitions ~ CMDsm", type = "ignition", climVar = "CMDsm",
               family = NULL, digestOfData = list(), modelAlgorithm = "randomForest"),
    "Other modelAlgorithms not implemented", fixed = TRUE)
})

test_that("runGLM.NB returns the ROC curve of a negative binomial fit", {
  ## ignitions rise with x; x = 12 has none although x = 11 has one. With fitted values increasing
  ## in x, the 9 presences and 11 absences make 99 pairs, and the only one ranked wrongly is
  ## (x = 11, x = 12): AUC = 98 / 99.
  dat <- data.frame(ignitions = c(rep(0, 10), 1, 0, 1, 2, 1, 3, 2, 4, 3, 5), x = 1:20)
  r <- suppressMessages(suppressWarnings(runGLM.NB(dat)))
  expect_equal(as.numeric(r$auc), 98 / 99)
  expect_identical(as.numeric(r$response), as.numeric(dat$ignitions > 0))
  expect_true(all(diff(as.numeric(r$predictor)) > 0))                    # fitted values rise with x
})
