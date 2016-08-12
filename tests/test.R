library(SpaDES)

## Not piecewise
  ## Poisson
    mySim <- simInit(
      times = list(start = 1, end = 1, timeunit = "year"),
      modules = list("fireSense_FrequencyFit"),
      paths = list(modulePath = " # replace with empty string instead"),
      params = list(fireSense_FrequencyFit = list(
        formula = formula(NB_FIRES_L ~ MDC_JUL * HW),
        family = "poisson",
        trace = 5,
        data = "dataFireSense_FrequencyFit")
      ),
      inputs = data.frame(
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
    #     formula = formula(NB_FIRES_L ~ MDC_JUL * HW),
    #     trace = 5,
    #     data = "dataFireSense_FrequencyFit")
    #   ),
    #   inputs = data.frame(
    #     files = "Z:/dataFireSense_FrequencyFit.rds",
    #     objectName = "dataFireSense_FrequencyFit",
    #     functions = "readRDS",
    #     package = "base",
    #     stringsAsFactors = FALSE)
    #   )

# Piecewise
  ## Poisson
    # mySim <- simInit(
    #   times = list(start = 1, end = 1, timeunit = "year"),
    #   modules = list("fireSense_FrequencyFit"),
    #   paths = list(modulePath = " # replace with empty string instead"),
    #   params = list(fireSense_FrequencyFit = list(
    #     formula = formula(NB_FIRES_L ~ MDC_JUL:HW + MDC_JUL:CN + MDC_JUL:DIST + MDC_JUL:O +
    #                         pw(MDC_JUL, K1):HW + pw(MDC_JUL, K2):CN + pw(MDC_JUL, K3):DIST + pw(MDC_JUL, K4):O - 1),
    #     family = "poisson",
    #     trace = 5,
    #     data = "dataFireSense_FrequencyFit")
    #   ),
    #   inputs = data.frame(
    #     files = "Z:/dataFireSense_FrequencyFit.rds",
    #     objectName = "dataFireSense_FrequencyFit",
    #     functions = "readRDS",
    #     package = "base",
    #     stringsAsFactors = FALSE)
    # )
  
  ## Negative binomial
    # mySim <- simInit(
    #   times = list(start = 1, end = 1, timeunit = "year"),
    #   modules = list("fireSense_FrequencyFit"),
    #   paths = list(modulePath = " # replace with empty string instead"),
    #   params = list(fireSense_FrequencyFit = list(
    #     formula = formula(NB_FIRES_L ~ MDC_JUL:HW + MDC_JUL:CN + MDC_JUL:DIST + MDC_JUL:O +
    #                       pw(MDC_JUL, K1):HW + pw(MDC_JUL, K2):CN + pw(MDC_JUL, K3):DIST + pw(MDC_JUL, K4):O - 1),
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
        data = "dataFireSense_FrequencyFit")
        objectName = "dataFireSense_FrequencyFit",

spades(mySim)
