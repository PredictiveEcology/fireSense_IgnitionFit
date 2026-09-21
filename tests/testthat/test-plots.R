## setupPlots() only needs objects with a predict() method, so exact linear models stand in for the
## xgboost ones and every plotted value can be worked out by hand.

plotInputs <- function() {
  dat <- scaleCovariates(makeIgnitionCovariates())
  ## y = 1 + 2 * CMDsm + 3 * youngAge exactly (fold 1), and twice that (fold 2)
  y <- 1 + 2 * dat$CMDsm + 3 * dat$youngAge
  covs <- as.data.frame(dat)[, c("CMDsm", "youngAge", "Pice_mar", "lightning")]
  list(dat = dat,
       models = list(Fold1 = stats::lm(y ~ ., data = cbind(y = y, covs)),
                     Fold2 = stats::lm(y ~ ., data = cbind(y = 2 * y, covs))))
}

test_that("setupPlots predicts along each covariate with the others at their mean", {
  inp <- plotInputs()
  plots <- setupPlots(inp$models, dat = inp$dat, igOrEsc = "ignition")
  expect_identical(names(plots), c("Climate", "Fuel"))

  clim <- as.data.frame(plots$Climate$data)
  fuel <- as.data.frame(plots$Fuel$data)
  ## CMDsm and lightning are the climate covariates; pixelID, year and ignitions are not plotted
  expect_setequal(unique(clim$variable), c("CMDsm", "lightning"))
  expect_setequal(unique(fuel$variable), c("youngAge", "Pice_mar"))

  ## youngAge is 0/1 with mean 0.5, so standardised it is +-0.9987: the x values run from
  ## floor(-9.987) / 10 = -1 to ceiling(9.987) / 10 = 1 in steps of 0.1, once per fold
  ya <- fuel[fuel$variable == "youngAge", ]
  expect_equal(sort(unique(ya$value)), seq(-1, 1, by = 0.1))
  expect_identical(nrow(ya), 2L * 21L)
  ## CMDsm spans 0..10, mean 5, sd 2.913: +-1.716 -> -1.8 .. 1.8
  cm <- clim[clim$variable == "CMDsm", ]
  expect_equal(range(cm$value), c(-1.8, 1.8))

  ## predictions: along youngAge, 1 + 3 * value (fold 1) and 2 + 6 * value (fold 2)
  expect_equal(sort(ya$predictedProb), sort(c(1 + 3 * seq(-1, 1, by = 0.1), 2 + 6 * seq(-1, 1, by = 0.1))))
  ## along CMDsm, 1 + 2 * value and 2 + 4 * value; along Pice_mar, which has no effect, 1 and 2
  expect_equal(range(cm$predictedProb), c(2 + 4 * -1.8, 2 + 4 * 1.8))
  expect_equal(sort(unique(round(fuel$predictedProb[fuel$variable == "Pice_mar"], 8))), c(1, 2))
  expect_identical(plots$Fuel$labels$y, "Predicted probability of ignition")
})

test_that("setupPlots limits the x axis to +-2 standard deviations", {
  inp <- plotInputs()
  data.table::set(inp$dat, 1L, "Pice_mar", 7.34)                 # an outlier, 7 sd from the mean
  plots <- setupPlots(inp$models["Fold1"], dat = inp$dat, igOrEsc = "escape")
  fuel <- as.data.frame(plots$Fuel$data)
  pm <- fuel[fuel$variable == "Pice_mar", ]
  ## Pice_mar is 0..96, mean 47.9, sd 28.0: the minimum is -1.71, and floor(-17.1) / 10 = -1.8
  expect_equal(range(pm$value), c(-1.8, 2))
  ## the sequence runs to ceiling(73.4) = 74 tenths; 20, 21, ..., 74 are all set to 2: 55 rows
  expect_identical(sum(pm$value == 2), 55L)
  expect_identical(plots$Climate$labels$y, "Predicted probability of escape")
})

test_that("plotPredictions plots one group of covariates in the colours given", {
  df <- data.table::data.table(
    value = c(-1, 0, 1, -1, 0, 1, -1, 0, 1),
    predictedProb = c(0.1, 0.2, 0.3, 0.5, 0.5, 0.5, 0.9, 0.8, 0.7),
    varFac = factor(rep(c("CMDsm", "Pice_mar", "youngAge"), each = 3),
                    levels = c("CMDsm", "Pice_mar", "youngAge"))
  )
  labels <- c(CMDsm = "Climate", Pice_mar = "Fuel", youngAge = "Fuel")
  colors <- c(CMDsm = "#FF0000", Pice_mar = "#00FF00", youngAge = "#0000FF")

  p <- plotPredictions(df, fuelOrClimate = "Fuel", labels = labels, jitter = 0.05, colors = colors,
                       fuelOrClimateInd = 2:3, igOrEsc = "escape")
  expect_identical(as.character(unique(p$data$varFac)), c("Pice_mar", "youngAge"))
  expect_identical(nrow(p$data), 6L)
  expect_identical(p$labels$x, "Scaled, centred")
  expect_identical(p$labels$y, "Predicted probability of escape")
  expect_identical(unname(vapply(p$layers, function(l) class(l$geom)[1], character(1))),
                   c("GeomPoint", "GeomPoint", "GeomSmooth"))   # points, jittered points, smooth

  ## the points are drawn at the data values, each covariate in its own colour
  pts <- suppressWarnings(ggplot2::ggplot_build(p))$data[[1]]   # loess warns on 3 points
  expect_equal(pts$x, c(-1, 0, 1, -1, 0, 1))
  expect_equal(pts$y, c(0.5, 0.5, 0.5, 0.9, 0.8, 0.7))
  expect_identical(pts$colour, rep(c("#00FF00", "#0000FF"), each = 3))

  clim <- plotPredictions(df, fuelOrClimate = "Climate", labels = labels, jitter = 0.05,
                          colors = colors, fuelOrClimateInd = 1L, igOrEsc = "ignition")
  expect_identical(unique(suppressWarnings(ggplot2::ggplot_build(clim))$data[[1]]$colour), "#FF0000")
  expect_equal(clim$data$predictedProb, c(0.1, 0.2, 0.3))
})
