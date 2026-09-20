## `dig` reaches runXGBOOST() as the cache key of the covariates. Before the fix, the module
## asked prepareCovariatesOuter() for `digestOfData$fireSense_ignitionCovariates`, which is not
## an element of what that function returns, so NULL was passed and the data took no part in the
## cache key: two different datasets shared a cache entry. dataDigest() takes the element that
## is actually there, `unscaledData`, and refuses to pass nothing.

test_that("dataDigest returns the digest prepareCovariatesOuter actually produces", {
  ## the shape of digestOfData: a named list with one element, `unscaledData`
  dig <- reproducible::.robustDigest(list(unscaledData = makeIgnitionCovariates()))
  expect_identical(names(dig), "unscaledData")

  expect_identical(dataDigest(list(digestOfData = dig)), dig[["unscaledData"]])
  ## the old, wrong key is not an element of it, so the bug is not silently reachable
  expect_null(dig[["fireSense_ignitionCovariates"]])
})

test_that("dataDigest distinguishes datasets that differ, and is stable for identical ones", {
  d1 <- makeIgnitionCovariates()
  d2 <- data.table::copy(d1)
  data.table::set(d2, 1L, "CMDsm", d2$CMDsm[1] + 1)

  key <- function(d) dataDigest(list(digestOfData = reproducible::.robustDigest(list(unscaledData = d))))
  expect_identical(key(d1), key(data.table::copy(d1)))
  expect_false(identical(key(d1), key(d2)))
})

test_that("dataDigest stops rather than cache the fit without the data in its key", {
  msg <- "no digest of the covariates"
  expect_error(dataDigest(list(digestOfData = NULL)), msg)
  expect_error(dataDigest(list(digestOfData = list())), msg)
  ## this is precisely the old behaviour: the wrong name, therefore NULL
  expect_error(dataDigest(list(digestOfData = list(somethingElse = "abc"))), msg)
  expect_error(dataDigest(list()), msg)
})

test_that("buildModel reaches runXGBOOST with the digest, not NULL", {
  ## buildModel() has no return value that shows `dig`, so it is checked where it lands: with a
  ## digest that does not match a cached fit, the fit is recomputed rather than restored.
  withr::local_options(list(reproducible.cachePath = withr::local_tempdir()))
  fit <- function(d, dig) {
    res <- NULL
    utils::capture.output(suppressWarnings(suppressMessages(
      res <- buildModel(covariates = data.table::copy(d), type = "ignition",
                        dig = dig, modelAlgorithm = "xgboost", nFolds = 5))))
    res
  }
  d1 <- scaleCovariates(makeIgnitionCovariates())
  d2 <- data.table::copy(d1)
  data.table::set(d2, NULL, "ignitions", rev(d1$ignitions))
  digOf <- function(d) dataDigest(list(digestOfData = reproducible::.robustDigest(list(unscaledData = d))))

  set.seed(42); f1 <- fit(d1, digOf(d1))
  set.seed(42); f2 <- fit(d2, digOf(d2))
  expect_false(isTRUE(all.equal(unname(aucPerFold(f1$rocs)), unname(aucPerFold(f2$rocs)))))
})

test_that("two different datasets do not share a cached xgboost fit", {
  ## the consequence of the bug: with dig = NULL the cache key of the fold models carried no
  ## information about the data, so a second, different dataset reused the first fit.
  cachePath <- withr::local_tempdir()
  withr::local_options(list(reproducible.cachePath = cachePath))
  fit <- function(d, dig) {
    res <- NULL
    utils::capture.output(suppressWarnings(suppressMessages(
      res <- runXGBOOST(data.table::copy(d), dig = dig, type = "ignition", nFolds = 5))))
    res
  }
  d1 <- scaleCovariates(makeIgnitionCovariates())
  d2 <- data.table::copy(d1)
  data.table::set(d2, NULL, "ignitions", rev(d1$ignitions))

  set.seed(42); f1 <- fit(d1, dataDigest(list(digestOfData = reproducible::.robustDigest(list(unscaledData = d1)))))
  set.seed(42); f2 <- fit(d2, dataDigest(list(digestOfData = reproducible::.robustDigest(list(unscaledData = d2)))))
  ## the second fit is not the first one, recovered from the cache
  expect_false(isTRUE(all.equal(unname(aucPerFold(f1$rocs)), unname(aucPerFold(f2$rocs)))))
})
