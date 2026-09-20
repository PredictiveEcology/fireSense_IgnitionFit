## The module run as a module: simInit() + spades() on the synthetic inputs. This is the only way
## to reach doEvent, Init, frequencyFitRun, buildModelsFitModels and .inputObjects, which need a
## simList.

runModule <- function(params = list(), objects = list(), end = 1, run = TRUE, drop = character()) {
  defaults <- list(
    fireSense_ignitionCovariates = makeIgnitionCovariates(),
    fireSense_escapeCovariates   = makeIgnitionCovariates(escapes = TRUE),
    ignitionFitRTM               = makeIgnitionFitRTM()
  )
  defaults[names(objects)] <- objects
  objects <- defaults
  objects <- objects[setdiff(names(objects), drop)]
  root <- withr::local_tempdir(.local_envir = parent.frame())
  paths <- testPaths
  paths$cachePath <- file.path(root, "cache")
  paths$outputPath <- file.path(root, "outputs")
  sim <- NULL
  utils::capture.output(suppressWarnings(suppressMessages({
    sim <- SpaDES.core::simInit(times = list(start = 1, end = end), modules = moduleName,
                                paths = paths, objects = objects,
                                params = stats::setNames(list(params), moduleName))
    if (run) sim <- SpaDES.core::spades(sim)
  })))
  sim
}

moduleEvents <- function(dt) {
  dt <- as.data.frame(dt)
  dt[dt$moduleName == moduleName & dt$eventType != ".inputObjects", c("eventTime", "eventType")]
}

test_that("init schedules checkData then run, and both processes are fitted and plotted", {
  sim <- runModule(params = list(.plots = "png"))

  done <- moduleEvents(SpaDES.core::completed(sim))
  expect_identical(done$eventType, c("init", "checkData", "run"))
  expect_equal(done$eventTime, c(1, 1, 1))
  expect_identical(nrow(moduleEvents(SpaDES.core::events(sim))), 0L)   # .runInterval = NA: nothing left

  for (nm in c("fireSense_IgnitionFitted", "fireSense_EscapeFitted")) {
    fitted <- sim[[nm]]
    expect_identical(names(fitted), c("modelList", "scaleData"))
    ml <- fitted$modelList
    expect_identical(names(ml), c("model", "rescales", "fittingRes", "lambdaRescaleFactor"))
    expect_identical(names(ml$model), c(paste0("Fold", 1:5), "rocs"))
    expect_equal(ml$fittingRes, 250)                       # cell size of ignitionFitRTM, m
    expect_equal(ml$lambdaRescaleFactor, 0.4)              # 400 rows / 1000 non-NA cells
    expect_null(ml$rescales)
    ## the known signal is recovered (see test-runXGBOOST.R)
    expect_equal(unname(aucPerFold(ml$model$rocs)), rep(1, 5))
    expect_setequal(stats::variable.names(ml$model$Fold1), c("CMDsm", "lightning", "Pice_mar", "youngAge"))
  }
  expect_identical(class(sim$fireSense_IgnitionFitted$modelList), "fireSense_IgnitionFit")
  expect_identical(class(sim$fireSense_EscapeFitted$modelList), "fireSense_EscapeFit")

  ## scaleData is what a user needs to standardise new data: the mean and sd of each covariate.
  ## mean(1:100) = 50.5; mean(2001:2004) = 2002.5; youngAge is half 0, half 1;
  ## sd(year) = sqrt(100 * (1.5^2 + 0.5^2 + 0.5^2 + 1.5^2) / 399) = 1.119434
  sc <- sim$fireSense_IgnitionFitted$scaleData
  expect_identical(sc$dimnames[[2]], c("pixelID", "year", "CMDsm", "youngAge", "Pice_mar", "lightning"))
  expect_equal(sc$`scaled:center`[c("pixelID", "year", "youngAge")],
               c(pixelID = 50.5, year = 2002.5, youngAge = 0.5))
  expect_equal(sc$`scaled:scale`[["year"]], sqrt(500 / 399))
  raw <- makeIgnitionCovariates()
  expect_equal(sc$`scaled:center`[["CMDsm"]], mean(raw$CMDsm))
  expect_equal(sc$`scaled:scale`[["Pice_mar"]], stats::sd(raw$Pice_mar))
  ## the escape table has the same covariates: `escapes` is not standardised
  expect_identical(sim$fireSense_EscapeFitted$scaleData$dimnames[[2]], sc$dimnames[[2]])

  ## new data standardised with scaleData predict the observed counts
  covs <- c("CMDsm", "lightning", "Pice_mar", "youngAge")
  newdata <- as.data.frame(raw)[, covs]
  for (col in covs) newdata[[col]] <- (newdata[[col]] - sc$`scaled:center`[[col]]) / sc$`scaled:scale`[[col]]
  pred <- predict(sim$fireSense_IgnitionFitted$modelList$model$Fold1, newdata)
  expect_equal(as.numeric(tapply(pred, raw$ignitions, mean)), c(0, 1, 2), tolerance = 0.1)

  ## the inputs are not changed by the fit
  expect_equal(as.data.frame(sim$fireSense_ignitionCovariates), as.data.frame(raw))
  expect_equal(as.data.frame(sim$fireSense_escapeCovariates),
               as.data.frame(makeIgnitionCovariates(escapes = TRUE)))

  ## one png per process
  figs <- list.files(SpaDES.core::figurePath(sim), recursive = TRUE, pattern = "[.]png$", full.names = TRUE)
  expect_length(figs, 2L)
  expect_length(grep("FuelClimate_Lightning_predicted_ignition_time-ordered_", basename(figs)), 1L)
  expect_length(grep("FuelClimate__predicted_escape_time-ordered_", basename(figs)), 1L)
  expect_true(all(file.size(figs) > 1000))
})

test_that("only the requested process is fitted, nothing is plotted without .plots, and run repeats", {
  sim <- runModule(params = list(whichProcessesToFit = "ignition", .plots = NA, .runInterval = 1),
                   end = 2, drop = "fireSense_escapeCovariates")
  expect_identical(names(sim$fireSense_IgnitionFitted$modelList$model), c(paste0("Fold", 1:5), "rocs"))
  expect_null(sim$fireSense_EscapeFitted)
  expect_length(list.files(SpaDES.core::figurePath(sim), recursive = TRUE), 0L)

  done <- moduleEvents(SpaDES.core::completed(sim))
  expect_identical(done$eventType, c("init", "checkData", "run", "run"))
  expect_equal(done$eventTime, c(1, 1, 1, 2))
  left <- moduleEvents(SpaDES.core::events(sim))
  expect_identical(left$eventType, "run")                  # the next one, past the end of the run
  expect_equal(left$eventTime, 3)
})

test_that(".runInitialTime sets when checkData and run happen", {
  sim <- runModule(params = list(.runInitialTime = 2), end = 2, run = FALSE)
  expect_identical(moduleEvents(SpaDES.core::events(sim))$eventType, "init")
  utils::capture.output(suppressWarnings(suppressMessages(
    sim <- SpaDES.core::spades(sim, events = list(fireSense_IgnitionFit = "init")))))
  sched <- moduleEvents(SpaDES.core::events(sim))
  expect_identical(sched$eventType, c("checkData", "run"))
  expect_equal(sched$eventTime, c(2, 2))
  expect_null(sim$fireSense_IgnitionFitted)
})

test_that("checkData stops when ignitionFitRTM has no nonNAs attribute", {
  rtm <- makeIgnitionFitRTM()
  attr(rtm, "nonNAs") <- NULL
  expect_error(runModule(objects = list(ignitionFitRTM = rtm), params = list(.plots = NA)),
               "nonNAs must be a non-empty/non-NULL numeric", fixed = TRUE)
  attr(rtm, "nonNAs") <- numeric(0)
  expect_error(runModule(objects = list(ignitionFitRTM = rtm), params = list(.plots = NA)),
               "nonNAs must be a non-empty/non-NULL numeric", fixed = TRUE)
})

test_that("checkData stops when whichProcessesToFit names neither process", {
  expect_error(runModule(params = list(whichProcessesToFit = "Ignition", .plots = NA)),
               "please review P(sim)$whichProcesesToFit", fixed = TRUE)
})

test_that(".inputObjects needs the covariates", {
  expect_error(runModule(drop = "fireSense_ignitionCovariates", run = FALSE),
               "this module does not produce data", fixed = TRUE)
})

test_that("rescaleVars = FALSE fits on the covariates as supplied and returns no scaleData", {
  sim <- runModule(params = list(whichProcessesToFit = "ignition", .plots = NA, rescaleVars = FALSE),
                   drop = "fireSense_escapeCovariates")
  expect_null(sim$fireSense_IgnitionFitted$scaleData)
  expect_equal(unname(aucPerFold(sim$fireSense_IgnitionFitted$modelList$model$rocs)), rep(1, 5))
  ## so the model predicts from unstandardised values
  raw <- makeIgnitionCovariates()
  pred <- predict(sim$fireSense_IgnitionFitted$modelList$model$Fold1,
                  as.data.frame(raw)[, c("CMDsm", "lightning", "Pice_mar", "youngAge")])
  expect_equal(as.numeric(tapply(pred, raw$ignitions, mean)), c(0, 1, 2), tolerance = 0.1)
})
