library(SpaDES)

## TODO: Expliquer que les coefficients du trace sont scalés et également que les bounds doivent être scalées (ou modif code).

# Define simulation parameters
times <- list(start = 1, end = 1, timeunit = "year")
modules <- list("fireSense_FrequencyFit")
paths <- list(
  modulePath = "~/Documents/GitHub/McIntire-lab/modulesPrivate/"
)


# Examples of model formula
formula <- n_fires ~ hw:MDC_07 + cn:MDC_07 + ot:MDC_07 + dt:MDC_07 +
  hw:pw(MDC_07, k_HW) + cn:pw(MDC_07, k_CN) + ot:pw(MDC_07, k_OT) + dt:pw(MDC_07, k_DT) - 1

# Define module parameters
parameters <- list(
  fireSense_FrequencyFit = list(
    formula = formula,
    family = quote(MASS::negative.binomial(theta = 1, link = "identity")),
    ub = list(coef = c(1, 1, 1, 1, 1, 1, 1, 1)),
    # start = c(1, 1),
    data = "dataFireSense_FrequencyFit",
    trace = 10, # Print progress every 10 iterations
    itermax = 10,
    nTrials = 2,
    nCores = 2
  )
)

# Define from where and how data will be loaded in the simList environment
inputs <- data.frame(
  objectName = "dataFireSense_FrequencyFit",
  file = "C:/Z/Contrats/Pessiere/DataInputs/dataFireSense_FrequencyFit_foudre.rds",
  fun = "readRDS",
  package = "base",
  loadTime = 1
)

# Create the simList
sim <- simInit(
  times = times,
  modules = modules,
  params = parameters,
  paths = paths,
  inputs = inputs
)

sim <- spades(sim)
sim$fireSense_FrequencyFitted

