## The module's metadata is its public contract: a project using this module binds
## to these object names and classes. Renaming or retyping one breaks every caller,
## which is exactly the class of change the raster -> terra migration makes, so it is
## worth asserting here rather than discovering downstream.
##
## When a change is deliberate, update this file in the same commit and bump the
## module version to match: removed, renamed or retyped is a MAJOR bump.

test_that("module metadata parses", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_type(md, "list")
  expect_identical(md$name, moduleName)
})

test_that("inputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  inputs <- stats::setNames(md$inputObjects$objectClass, md$inputObjects$objectName)
  expect_identical(
    inputs[order(names(inputs))],
    c(climateVariablesForFire      = "list",
      fireSense_ignitionCovariates = "data.frame",
      fireSense_ignitionFormula    = "character",
      ignitionFitRTM               = "SpatRaster")
  )
})

test_that("outputs are the expected names and classes", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  outputs <- stats::setNames(md$outputObjects$objectClass, md$outputObjects$objectName)
  expect_identical(
    outputs[order(names(outputs))],
    c(fireSense_EscapeFitted   = "fireSense_EscapeFit",
      fireSense_IgnitionFitted = "fireSense_IgnitionFit")
  )
})

test_that("parameters are the expected names", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  expect_identical(
    sort(md$parameters$paramName),
    sort(c(".plotInitialTime", ".plots", ".runInitialTime", ".runInterval",
           ".saveInitialTime", ".saveInterval", ".seed", ".studyAreaName", ".useCache",
           "crossValType", "escapeFamily", "ignitionFamily", "modelAlgorithm",
           "plot_fuelBiomassPerPrediction", "rescaleVars", "whichProcessesToFit"))
  )
})

test_that("parameters have the expected classes and defaults", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  p <- md$parameters
  classes <- stats::setNames(unlist(p$paramClass), p$paramName)
  expect_identical(
    classes[order(names(classes))],
    c(.plotInitialTime = "numeric", .plots = "character", .runInitialTime = "numeric",
      .runInterval = "numeric", .saveInitialTime = "numeric", .saveInterval = "numeric",
      .seed = "list", .studyAreaName = "character", .useCache = "logical",
      crossValType = "character", escapeFamily = "function, character",
      ignitionFamily = "function, character", modelAlgorithm = "character",
      plot_fuelBiomassPerPrediction = "numeric", rescaleVars = "logical",
      whichProcessesToFit = "character")
  )

  default <- function(name) p$default[[match(name, p$paramName)]]
  expect_identical(default("crossValType"), c("time-ordered", "crossValidation"))
  expect_identical(default("whichProcessesToFit"), c("ignition", "escape"))
  expect_identical(default("modelAlgorithm"), "xgboost")
  expect_identical(default("rescaleVars"), TRUE)
  expect_identical(default("escapeFamily"), quote(binomial(link = "logit")))
  expect_identical(default("ignitionFamily"), quote(poisson(link = "log")))
  expect_identical(default(".plots"), "screen")
  expect_identical(default(".useCache"), FALSE)
  for (nm in c("plot_fuelBiomassPerPrediction", ".plotInitialTime", ".seed"))
    expect_null(default(nm))
  for (nm in c(".runInterval", ".saveInitialTime", ".saveInterval", ".studyAreaName"))
    expect_true(is.na(default(nm)))
  ## the only parameter with bounds
  i <- match("plot_fuelBiomassPerPrediction", p$paramName)
  expect_equal(c(p$min[[i]], p$max[[i]]), c(1, 10))
})

test_that("every package the tests rely on is in reqdPkgs, so CI cannot skip them", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  pkgs <- sub("\\s*\\(.*$", "", sub("@.*$", "", sub("^.*/", "", unlist(md$reqdPkgs))))
  expect_true(all(c("xgboost", "SHAPforxgboost", "caret", "pROC", "MASS", "data.table", "terra",
                    "ggplot2", "ggpubr", "fireSenseUtils", "reproducible", "SpaDES.core") %in% pkgs))
})
