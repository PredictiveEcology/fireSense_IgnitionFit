## Loading mirai loads nanonext, which starts ~30 NNG threads the moment the namespace
## loads. SpaDES.core loads every module's reqdPkgs at simInit, so every job process then
## carries those threads, and any later fork (parallel::mcMap / mclapply, e.g. in
## fireSenseUtils::bufferToArea via fireSense_dataPrepFit) hangs forever on a lock the
## threads held. All mirai use in this module is commented out, so it must not be listed.
test_that("reqdPkgs does not list mirai", {
  exprs <- parse(testthat::test_path("..", "..", "fireSense_IgnitionFit.R"), keep.source = FALSE)
  dm <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("defineModule")), exprs)
  expect_length(dm, 1L)
  pkgs <- unlist(eval(dm[[1]][[3]]$reqdPkgs, baseenv()))
  pkgNames <- sub("\\s*\\(.*$", "", sub("@.*$", "", sub("^.*/", "", pkgs)))
  expect_false("mirai" %in% pkgNames)
})
