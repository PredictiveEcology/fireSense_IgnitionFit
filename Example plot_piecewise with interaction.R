library(dplyr)
library(ggplot2)
library(MASS)
library(parallel)
library(SpaDES)

modulePath <- "~/Documents/GitHub/McIntire-lab/modulesPrivate/"
start <- end <- 1

# Define simulation parameters
times <- list(start = start, end = end, timeunit = "year")
modules <- list("fireSense_FrequencyPredict")
paths <- list(
  modulePath = modulePath
)

model <- structure(
  list(
    formula = n_fires ~ hw:weather + cn:weather + hw:pw(weather, k_HW) + cn:pw(weather, k_CN) - 1, 
    family = negative.binomial(.42, link = "identity"),
    coef = structure(c(1e-23, 9e-05,  1e-03, 7e-03),
                     .Names = c("hw:weather", "cn:weather", "hw:pw(weather, k_HW)", "cn:pw(weather, k_CN)")),
    knots = structure(c(70, 33), .Names = c("k_HW", "k_CN")),
    theta = .42
  ), .Names = c("formula", "family", "coef", "knots", "theta"), class = "fireSense_FrequencyFit"
)

plotData <- data_frame(weather = rep(1:100, 2),
                       hw = rep(c(1, 0), each = 100),
                       cn = rep(c(0, 1), each = 100), 
                       group = rep(c("hw", "cn"), each = 100)
)

# Define module parameters
parameters <- list(
  fireSense_FrequencyPredict = list(
    modelName = "fireSense_FrequencyFitted",
    data = "plotData",
    modelName = "model"
  )
)

# Objects to pass from the global environment to the simList environment
objects <- c("plotData", "model")

# Create the simList
sim <- simInit(
  times = times, 
  params = parameters, 
  modules = modules, 
  objects = objects, 
  paths = paths
)

sim <- spades(sim)

# Prepare data
plotData <- bind_cols(plotData, list(predict = sim$fireSense_FrequencyPredicted[[as.character(start)]])) 

# Predicted number of fires as a function of two covariates
p <- ggplot(data = plotData) + theme_bw()
p <- p + geom_line(aes(x = weather, y = predict, group = group, color = group))
p
