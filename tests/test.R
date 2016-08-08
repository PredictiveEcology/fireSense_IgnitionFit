library(SpaDES)

mySim <- simInit(
  times = list(start = 1, end = 1, timeunit = "year"),
  modules = list("fireSense_FrequencyFit"),
  paths = list(modulePath = " # replace with empty string instead"),
  params = list(fireSense_FrequencyFit = list(
    #formula = formula(NB_FIRES_L ~ MDC_JUL * HW),
#    formula = formula(NB_FIRES_L ~ pw(MDC_JUL, K) * HW),
    formula = formula(NB_FIRES_L ~ MDC_JUL:HW + MDC_JUL:CN + MDC_JUL:DIST + MDC_JUL:O +
                     pw(MDC_JUL, BP1):HW + pw(MDC_JUL, BP2):CN + pw(MDC_JUL, BP3):DIST + pw(MDC_JUL, BP4):O - 1),
    family = "negative.binomial",
    trace = 5,
    data = "dataFireSense_FrequencyFit"
  )),
  inputs = data.frame(
    files = "Z:/dataFireSense_FrequencyFit.rds",
    objectName = "dataFireSense_FrequencyFit",
    functions = "readRDS",
    package = "base",
    stringsAsFactors = FALSE)
  )

spades(mySim)
