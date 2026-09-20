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
    c(fireSense_ignitionCovariates = "data.frame",
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
    sort(c(".plots", ".runInitialTime", ".runInterval", ".seed", ".useCache",
           "crossValType", "modelAlgorithm", "rescaleVars", "whichProcessesToFit"))
  )
})

test_that("parameters have the expected classes and defaults", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  p <- md$parameters
  classes <- stats::setNames(unlist(p$paramClass), p$paramName)
  expect_identical(
    classes[order(names(classes))],
    c(.plots = "character", .runInitialTime = "numeric", .runInterval = "numeric",
      .seed = "list", .useCache = "logical", crossValType = "character",
      modelAlgorithm = "character", rescaleVars = "logical",
      whichProcessesToFit = "character")
  )

  default <- function(name) p$default[[match(name, p$paramName)]]
  expect_identical(default("crossValType"), c("time-ordered", "crossValidation"))
  expect_identical(default("whichProcessesToFit"), c("ignition", "escape"))
  expect_identical(default("modelAlgorithm"), "xgboost")
  expect_identical(default("rescaleVars"), TRUE)
  expect_identical(default(".plots"), "screen")
  expect_identical(default(".useCache"), FALSE)
  expect_null(default(".seed"))
  expect_true(is.na(default(".runInterval")))
  ## no parameter has bounds any more
  expect_true(all(vapply(p$min, function(x) is.na(x) || is.null(x), logical(1))))
})

test_that("every package the tests rely on is in reqdPkgs, so CI cannot skip them", {
  md <- SpaDES.core::moduleMetadata(module = moduleName, path = modulePath)
  pkgs <- sub("\\s*\\(.*$", "", sub("@.*$", "", sub("^.*/", "", unlist(md$reqdPkgs))))
  expect_true(all(c("xgboost", "SHAPforxgboost", "caret", "pROC", "data.table", "terra",
                    "ggplot2", "ggpubr", "fireSenseUtils", "reproducible", "SpaDES.core") %in% pkgs))
})
