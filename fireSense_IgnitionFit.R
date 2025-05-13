defineModule(sim, list(
  name = "fireSense_IgnitionFit",
  description = paste("Fit statistical models that can be used to parameterize (calibrate)",
                      "the fire ignition component of landscape fire models (e.g. fireSense)."),
  keywords = c("fire frequency", "optimization", "additive property", "poisson",
               "negative binomial", "fireSense"),
  authors = c(
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut")),
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = "aut"),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(),
  version = list(fireSense_IgnitionFit = "1.0.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_IgnitionFit.Rmd"),
  loadOrder = list(after = "fireSense_dataPrepFit"),
  reqdPkgs = list("data.table", "dplyr",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.0.5.9090)",
                  "glmmTMB",
                  "ggplot2", "ggpubr", "MASS", "magrittr",
                  "numDeriv", "parallel", "parallelly",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/reproducible@AI (>= 2.1.2)",
                  #TODO correct this when reproducible is merged - it is due to cache(predict)
                  "RhpcBLASctl",
                  "PredictiveEcology/SpaDES.core@development (>= 2.0.2.9006)", "terra"),
  parameters = bindrows(
    defineParameter("escapeFamily", "function, character", default = quote(binomial(link = "logit")),
                    desc = paste("a family function (must be wrapped with `quote()`) or a",
                                 "character string naming a family function.",
                                 "Only the negative binomial has been implemented",
                                 "For additional details see `?family`. This was formerly ",
                                 "`quote(MASS::negative.binomial(theta = 1, link = 'identity'))`.")),
    defineParameter("ignitionFamily", "function, character", default = quote(poisson(link = "log")),
                    desc = paste("a family function (must be wrapped with `quote()`) or a",
                                 "character string naming a family function.",
                                 "Only the negative binomial has been implemented",
                                 "For additional details see `?family`. This was formerly ",
                                 "`quote(MASS::negative.binomial(theta = 1, link = 'identity'))`.")),
    defineParameter("plot_fuelBiomassPerPrediction", "numeric", NULL, 1, 10,
                    desc = paste("when generating plots of climate x fuel class, the log of biomass (g/m2) for which",
                                 "to generate predictions across a gradient of climate values.",
                                 "If supplied, it will override any values in `sim$ignitionFitRTM$meanForestB`")),
    defineParameter("rescaleVars", "logical", default = TRUE,
                    desc = paste("Attempt to rescale variables? If `rescalers` is defined,",
                                 "use it to rescale variables as `var / rescalers['var']`. ",
                                 "Otherwise, `scale()` will be used to rescale variables to `[0,1]`,",
                                 "if they are not already within this range.")),
    defineParameter("whichProcessesToFit", "character", c("ignition", "escape"), NA, NA,
                    "which processes to fit: ignition, escape, or both (the default)"),
    defineParameter(".plots", "character", default = "screen",
                    desc = "See ?Plots. There are a few plots that are made within this module, if set."),
    defineParameter(".plotInitialTime", "numeric", default = NULL,
                    desc = "when to do plot"),
    defineParameter(".runInitialTime", "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start time of the simulation."),
    defineParameter(".runInterval", "numeric", default = NA,
                    desc = paste("optional. Interval between two runs of this module,",
                                 "expressed in units of simulation time.",
                                 "By default, NA, which means that this module only runs once per simulation.")),
    defineParameter(".saveInitialTime", "numeric", default = NA,
                    desc = "optional. When to start saving output to a file."),
    defineParameter(".saveInterval", "numeric", default = NA,
                    desc = "optional. Interval between save events."),
    defineParameter(".seed", "list", NULL, NA, NA,
                    paste("Named list of seeds to use for each event (names).",
                          "E.g., `list('init' = 123)` will `set.seed(123)`",
                          "at the start of the init event and unset it at the end.",
                          "Defaults to `NULL`, meaning that no seeds will be set.")),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used. If NA, a hash of `studyAreaLarge` will be used."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules,",
                                 "where stochasticity and time are not relevant."))
  ),
  inputObjects = bindrows(
    expectsInput("climateVariablesForFire", "list",
                 desc = paste("The column name in the `fireSense_ignitionCovariates` that is climate,",
                              "in a named list, .e.g. `climateVariablesForFire = list('ignition' = 'MDC')`")),
    expectsInput("fireSense_ignitionCovariates", "data.frame",
                 desc = "table of aggregated ignition covariates with annual ignitions"),
    expectsInput("flammableRTM", "SpatRaster", sourceURL = NA,
                 "RTM without ice/rocks/urban/water. Flammable map with 0 and 1."),
    expectsInput("ignitionFitRTM", "SpatRaster",
                 desc = paste("A (template) raster with information with regards to the spatial",
                              "resolution and geographical extent of `fireSense_ignitionCovariates`.",
                              "Used to pass this information onto `fireSense_ignitionFitted`",
                              "Needs to have number of non-NA cells as attribute:",
                              "(`ignitionFitRTM@data@attributes$nonNAs`), and optionally,",
                              "`ignitionFitRTM@data@attributes$meanForestB`")),
    expectsInput("fireSense_ignitionFormula", "character",
                 desc = "formula - as a character - describing the model to be fitted."),
  ),
  outputObjects = bindrows(
    createsOutput("fireSense_IgnitionFitted", "fireSense_IgnitionFit",
                  desc = "A fitted model object of class `fireSense_IgnitionFit`.")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_IgnitionFit = function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {
      #the reason this is not scheduling Init is because it is only performing sanity checks on data
      # which is created during other modules non-init events. And since Inits are all scheduled first...
      sim <- scheduleEvent(sim, P(sim)$.runInitialTime, moduleName, "checkData", eventPriority = 2)

      sim <- scheduleEvent(sim, P(sim)$.runInitialTime, moduleName, "run", eventPriority = 5.11)

      if (!is.na(P(sim)$.saveInitialTime)) {
        sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, moduleName, "save", .last())
      }
    },
    checkData = {

      sim <- Init(sim)
    },
    run = {
      sim <- frequencyFitRun(sim)

      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "run")
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )

  return(invisible(sim))
}

### template initialization
Init <- function(sim) {

  if ("ignition" %in% P(sim)$whichProcessesToFit) {
    if (is.null(attributes(sim$ignitionFitRTM)$nonNAs) ||
        length(attributes(sim$ignitionFitRTM)$nonNAs) == 0) {
      stop("sim$ignitionFitRTM@data@attributes$nonNAs must be a non-empty/non-NULL numeric")
    }
    checkData(formula = sim$fireSense_ignitionFormula,
              covariates = sim$fireSense_ignitionCovariates, type = "ignition")
  }
  if ("escape" %in% P(sim)$whichProcessesToFit) {
    checkData(sim$fireSense_escapeFormula, sim$fireSense_escapeCovariates, "escape")
  }
  if (!any(c("ignition", "escape") %in% P(sim)$whichProcessesToFit)) {
    stop("please review P(sim)$whichProcesesToFit...ensure lower-case")
  }

  return(invisible(sim))
}

frequencyFitRun <- function(sim) {

  if ("ignition" %in% P(sim)$whichProcessesToFit) {


    ignitionData <- prepareCovariates(formula = sim$fireSense_ignitionFormula,
                                      covariates = sim$fireSense_ignitionCovariates,
                                      rescaleVars =P(sim)$rescaleVars)
    #use only a sample of zeroes...
    # zeroes <- ignitionData$covariates[ignitions == 0]
    # nonzeroes <- ignitionData$covariates[ignitions > 0]
    # #take 10 times more zeroes than igs
    # sampleZeroes <- zeroes[sample(nrow(zeroes), size = nrow(nonzeroes) * 10, replace = FALSE)]
    # igSample <- rbind(nonzeroes, sampleZeroes)

    ignitionModel <- buildModel(covariates = ignitionData$covariates,
                                climVar = sim$climateVariablesForFire$ignition,
                                formula= ignitionData$formula, type = "ignition",
                                family = P(sim)$ignitionFamily)
    #ignition specific
    origNoPix <- attributes(sim$ignitionFitRTM)$nonNAs   ## nrow(preSampleData) in eg above
    finalNoPix <- nrow(ignitionData$fireSense_ignitionCovariates)     ## nrow(postSampleData) in eg above
    lambdaRescaleFactor <- finalNoPix/origNoPix

    modelList <- list(
      model = ignitionModel$bestModel,
      rescales = ignitionData$ignitionRescalers,
      fittingRes = res(sim$ignitionFitRTM)[1],
      lambdaRescaleFactor = lambdaRescaleFactor)
    sim$fireSense_IgnitionFitted <- modelList
    class(sim$fireSense_IgnitionFitted) <- "fireSense_IgnitionFit"

    if (anyPlotting(P(sim)$.plots)) {
      IgEscapePlots(dt = ignitionData$covariates, bestModel = ignitionModel,
                    climVar = sim$climateVariablesForFire$ignition,
                    rescalers = ignitionData$ignitionRescalers,
                    fsProcess = "ignition", family =  P(sim)$ignitionFamily,
                    plotBiomass = P(sim)$plot_fuelBiomassPerPrediction,
                    ignitionFitRTM = sim$ignitionFitRTM,
                    studyAreaName = P(sim)$.studyAreaName,
                    oPath = outputPath(sim))
    }
  }

  if ("escape" %in% P(sim)$whichProcessesToFit) {
    escapeData <- prepareCovariates(formula = sim$fireSense_escapeFormula,
                                    covariates = sim$fireSense_escapeCovariates,
                                    rescaleVars = P(sim)$rescaleVars)

    escapeModel <- buildModel(covariates = escapeData$covariates,
                              climVar = sim$climateVariablesForFire$ignition,
                              formula= escapeData$formula, type = "escape",
                              family = P(sim)$escapeFamily)

    if (anyPlotting(P(sim)$.plots)) {

      IgEscapePlots(dt = escapeData$covariates, bestModel = escapeModel,
                    climVar = sim$climateVariablesForFire$ignition,
                    fsProcess = "escape", family =  P(sim)$escapeFamily,
                    rescalers = escapeData$ignitionRescalers,
                    plotBiomass = P(sim)$plot_fuelBiomassPerPrediction,
                    ignitionFitRTM = sim$ignitionFitRTM,
                    studyAreaName = P(sim)$.studyAreaName,
                    oPath = outputPath(sim))
    }
  }

  return(invisible(sim))
}


expit <- function(x) 1/(1 + exp(-x)) # inverse logit function; used below

messageFormulaFn <- function(form) {
  gsub(" {2,100}", " ", paste0(format(form), collapse = ""))
}

checkData <- function(formula, covariates, type) {
  if (!"pixelID" %in% colnames(covariates)) {
    stop("covariates for ", type, " must have a 'pixelID' column")
  }

  if (is.empty.model(as.formula(formula, env = .GlobalEnv))) {
    stop("formula for ", type, " describes an empty model.")
  }
}

prepareCovariates <- function(formula, covariates, rescaleVars) {

  formula <- as.formula(formula, env = .GlobalEnv)
  terms <- terms.formula(formula)

  covariates <- copy(setDT(covariates))

  if (attr(terms, "response")) {
    y <- formula[[2L]]
  } else {
    stop("Incomplete formula, the LHS is missing.")
  }

  if (any(c("year", "yr") %in% tolower(names(covariates)))) {
    xvar <- intersect(c("year", "yr"), tolower(names(covariates)))
  } else {
    xvar <- rows #TODO what is this?
  }

  if (rescaleVars) {
    # rescalers <- abs(sapply(covariates[, .SD, .SDcol = toRescale], FUN = max))
    message("Variables outside of [0,10] range will be rescaled to [0,10]")

    toRescale <- setdiff(names(covariates),
                         c("pixelID", "ignitions", "year", "yearChar"))
    rescalers <- sapply(covariates[, .SD, .SDcol = toRescale], max)
    needRescale <- sapply(rescalers, FUN = function(x) !inRange(x, 0, 10))
    cols <- names(rescalers)[which(needRescale)]
    message("rescaling the following variables: ", paste(cols, collapse = ", "))
    ignitionRescalers <- 10^floor(log10(abs(rescalers[cols])))
    covariates <- rescaleVarsByMagnitude(covariates, ignitionRescalers)
  } else {
    ignitionRescalers <- NULL #so that fire fireSense_IgnitionFit can add it
  }

  return(list(covariates = covariates,
              formula = formula,
              ignitionRescalers = ignitionRescalers,
              xvar = xvar))
}

buildModel <- function(covariates, formula,  type = "ignition",
                       climVar = sim$climateVariablesForFire$ignition,
                       family) {

  # convert to data.table --> easier to work with
  m <- as.data.table(covariates)

  # Run 3 models ... the full model, climate only,  and NULL model
  forms <- list()
  forms[["full"]] <- as.formula(formula, env = .GlobalEnv)
  # drop interactions
  formChar <- as.character(forms[["full"]])
  terms <- lapply(formChar, function(x) {
    allTermsNoMinus <- strsplit(x, " *\\- *")[[1]]
    allTermsNoMinus <- lapply(allTermsNoMinus, function(y) {
      allTerms <- strsplit(y, " *\\+ *")[[1]]
      noInteractions <- grep(":", allTerms, value = TRUE, invert = TRUE)
      paste0(noInteractions, collapse = " + ")
    })
    paste0(allTermsNoMinus, collapse = " - ")
  })

  forms[["NoInteractions"]] <- as.formula(paste0(terms[c(2,1,3)], collapse = " "), env = .GlobalEnv)
  forms[["climateOnly"]] <-  as.formula(paste0(
    paste0(terms[c(2,1,3)], collapse = " "),
    paste0(" + ", climVar)),
    env = .GlobalEnv)

  varToInt <- "ignitions"
  if (type == "escape") {
    varToInt <- c("ignitions", "escapes")
  }

  for (i in c("year", varToInt))
    set(m, NULL, i, as.integer(m[[i]]))
  for (i in c("yearChar"))
    set(m, NULL, i, factor(m[[i]]))

  system.time({
    mods <- Map(
      nam = names(forms), form = forms,
      MoreArgs = list(dat = m, family = family, type = type),
      f = function(form, nam, dat, family, type) {
        en <- new.env(parent = .GlobalEnv)
        if (type == "ignition") {
          ziform <- as.formula(paste0("~", paste0(climVar, collapse = "+")), env = en)
        } else {
          ziform <- as.formula(~0, env = en)
        }
        objNames <- c("dat", "family", "form", "ziform", "nam")
        objs <- mget(objNames)
        dig <- en$dig <- .robustDigest(objs)
        list2env(objs, envir = en)
        message("Running glmmTMB with Zero-Inflated, Mixed effect, Poisson, using:\n",
                messageFormulaFn(form))


        out <- local({
          glmmTMB(form, data = dat,
                  ziformula = ziform,
                  family = eval(family)) |>
            ## Use .cacheExtra: there are lots of arguments to glmmTMB that seemed to be "always different"
            #TODO: will nam be an issue if it is identical for escape and ignition models?
            Cache(.functionName = paste0("glmmTMB_for", "_", nam),
                  omitArgs = formalArgs(glmmTMB),
                  .cacheExtra = dig)},
          envir = en)

        # a <- identifyEnvs(out, en)
        out
      })
  })

  AICs <- sapply(mods, AIC)
  ## even if the AIC is <2 better, should take simpler model;
  ## in tests, turned many to non-significant when had interactions
  whBest <- which.min(c(AICs[["full"]] + 2, AICs[["NoInteractions"]], AICs[["climateOnly"]]))
  # whBest <- 1
  bestModel <- mods[[whBest]]
  messageColoured("Best model is:\n", messageFormulaFn(bestModel$call$formula), colour = "magenta")
  summ <- summary(bestModel)

  return(bestModel)
}


.inputObjects <- function(sim) {
  if (!suppliedElsewhere("fireSense_ignitionCovariates", sim)) {
    stop("this module does not produce data - consider usin the module 'PredictiveEcology/fireSense_dataPrepFit'")
  }

  if (!suppliedElsewhere("climateVariablesForFire", sim)) {
    #in most cases this will be supplied with data - this is to make some changes backwards compatible
    sim$climateVariablesForFire <- list("ignition" = "MDC")

    if (!"MDC" %in% names(sim$fireSense_ignitionCovariates)) {
      stop("please supply climateVariablesForFire")
    }
  }

  return(sim)
}
