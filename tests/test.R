library(SpaDES)

## Not piecewise
  ## Poisson
    # mySim <- simInit(
    #   times = list(start = 1, end = 1, timeunit = "year"),
    #   modules = list("fireSense_FrequencyFit"),
    #   paths = list(modulePath = " # replace with empty string instead"),
    #   params = list(fireSense_FrequencyFit = list(
    #     formula = formula(NB_FIRES_L ~ MDC_JUL - 1),
    #     family = poisson,
    #     trace = 5,
    #     data = "dataFireSense_Frequency")
    #   ),
    #   inputs = data.frame(
    #     files = "Z:/dataFireSense_Frequency.rds",
    #     objectName = NA,
    #     functions = "readRDS",
    #     package = "base",
    #     stringsAsFactors = FALSE)
    #   )

  ## Negative binomial
    # mySim <- simInit(
    #   times = list(start = 1, end = 1, timeunit = "year"),
    #   modules = list("fireSense_FrequencyFit"),
    #   paths = list(modulePath = " # replace with empty string instead"),
    #   params = list(fireSense_FrequencyFit = list(
    #     formula = formula(NB_FIRES_L ~ MDC_JUL * HW - 1),
    #     trace = 5,
    #     data = "dataFireSense_Frequency")
    #   ),
    #   inputs = data.frame(
    #     files = "Z:/dataFireSense_Frequency.rds",
    #     objectName = NA,
    #     functions = "readRDS",
    #     package = "base",
    #     stringsAsFactors = FALSE)
    #   )

## Piecewise
  ## Poisson
    mySim <- simInit(
      times = list(start = 1, end = 1, timeunit = "year"),
      modules = list("fireSense_FrequencyFit"),
      paths = list(modulePath = " # replace with empty string instead"),
      params = list(fireSense_FrequencyFit = list(
        formula = formula(NB_FIRES_L ~ MDC_JUL:HW + MDC_JUL:CN + MDC_JUL:D + MDC_JUL:O +
                            pw(MDC_JUL, K1):HW + pw(MDC_JUL, K2):CN + pw(MDC_JUL, K3):D + pw(MDC_JUL, K4):O - 1),
        family = "poisson",
        trace = 5,
        data = "dataFireSense_Frequency")
      ),
      inputs = data.frame(
        files = "Z:/dataFireSense_Frequency.rds",
        objectName = NA,
        functions = "readRDS",
        package = "base",
        stringsAsFactors = FALSE)
    )
    
  ## Negative binomial
    # mySim <- simInit(
    #   times = list(start = 1, end = 1, timeunit = "year"),
    #   modules = list("fireSense_FrequencyFit"),
    #   paths = list(modulePath = " # replace with empty string instead"),
    #   params = list(fireSense_FrequencyFit = list(
    #     formula = formula(NB_FIRES_L ~ MDC_JUL:HW + MDC_JUL:CN + MDC_JUL:D + MDC_JUL:O +
    #                       pw(MDC_JUL, K1):HW + pw(MDC_JUL, K2):CN + pw(MDC_JUL, K3):D + pw(MDC_JUL, K4):O - 1),
    #     family = MASS::negative.binomial(theta = 1, link = identity),
    #     ub = list(beta = rep(1e-1, 8)),
    #     trace = 5,
    #     data = "dataFireSense_Frequency")
    #   ),
    #   inputs = data.frame(
    #     files = "Z:/dataFireSense_Frequency.rds",
    #     objectName = NA,
    #     functions = "readRDS",
    #     package = "base",
    #     stringsAsFactors = FALSE)
    # )
#mySim$dataFireSense_Frequency$MDC_JUL <- -mySim$dataFireSense_Frequency$MDC_JUL ## Handling negative variables

spades(mySim)
