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
                  "glmmTMB", "mirai",
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
    createsOutput("fireSense_EscapeFitted", "fireSense_EscapeFit",
                  desc = "A fitted model object of class `fireSense_EscapeFit`"),
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
      sim <- scheduleEvent(sim, P(sim)$.runInitialTime, moduleName, "checkData")

      sim <- scheduleEvent(sim, P(sim)$.runInitialTime, moduleName, "run")

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

  #TODO: make cache smart by digest args in advance
  whichProcessesToFit <- P(sim)$whichProcessesToFit
  names(whichProcessesToFit) <- igOrEscNames(whichProcessesToFit, post = "Fitted", case = "sen")

  mir <- list()
  on.exit(lapply(mir, stop_mirai))
  fireSense_FittedModels <-
    Map(igOrEsc = whichProcessesToFit,
        function(igOrEsc) {
          formulaHere <- igOrEscNames(igOrEsc, post = "Formula") # paste0("fireSense_", igOrEsc, "Formula")
          covariatesHere <- igOrEscNames(igOrEsc, post = "Covariates") # paste0("fireSense_", igOrEsc, "Covariates")
          objsNeeded <- c(covariatesHere, formulaHere)
          modelFamily <- igOrEscNames(igOrEsc, pre = "", post = "Family") # paste0(igOrEsc, "Family")
          digestOfData <- .robustDigest(mget(objsNeeded, envir = envir(sim)))

          data <- prepareCovariates(formula = sim[[formulaHere]],
                                    covariates = sim[[covariatesHere]],
                                    rescaleVars =P(sim)$rescaleVars) |>
            Cache(omitArgs = c("covariates", "formula"), .cacheExtra = digestOfData)


          family <- P(sim)[[modelFamily]]
          modelHere <- buildModel(covariates = data$covariates,
                                  climVar = sim$climateVariablesForFire$ignition,
                                  formula= data$formula, type = "ignition",
                                  family = family) |>
            Cache(omitArgs = c("covariates", "formula"), .cacheExtra = c(digestOfData, 1))

          #ignition specific
          origNoPix <- attributes(sim$ignitionFitRTM)$nonNAs   ## nrow(preSampleData) in eg above
          finalNoPix <- nrow(data$covariates)     ## nrow(postSampleData) in eg above
          lambdaRescaleFactor <- finalNoPix/origNoPix

          modelList <- list(
            model = modelHere,
            rescales = data$ignitionRescalers,
            fittingRes = res(sim$ignitionFitRTM)[1],
            lambdaRescaleFactor = lambdaRescaleFactor)

          class(modelList) <- igOrEscNames(igOrEsc, post = "Fit", case = "Title")

          if (anyPlotting(P(sim)$.plots)) {
            library(mirai)
            try(daemons(2, dispatcher = FALSE), silent = TRUE) # this is ignored the 2nd time
            message("Plotting ", igOrEsc, "...")
            argsForDigest <- list(climVar = sim$climateVariablesForFire[[igOrEsc]],
                                  # rescalers = data$ignitionRescalers,
                                  fsProcess = igOrEsc, family = P(sim)[[modelFamily]],
                                  plotBiomass = P(sim)$plot_fuelBiomassPerPrediction,
                                  ignitionFitRTM = sim$ignitionFitRTM,
                                  studyAreaName = P(sim)$.studyAreaName,
                                  oPath = outputPath(sim), IgEscapePlots = IgEscapePlots,
                                  digestOfData = digestOfData)
            digestForPlot <- .robustDigest(argsForDigest)
            argsAll <- append(list(data = data, bestModel = modelHere, digestForPlot = digestForPlot),
                              argsForDigest)

            mir[[igOrEsc]] <<- mirai(.expr = {
              library(glmmTMB)
              library(reproducible)
              IgEscapePlots(dt = data$covariates, bestModel = bestModel,
                            climVar = climVar,
                            rescalers = data$ignitionRescalers,
                            fsProcess = fsProcess, family = family,
                            plotBiomass = plotBiomass,
                            ignitionFitRTM = ignitionFitRTM,
                            studyAreaName = studyAreaName,
                            oPath = oPath) |>
                Cache(omitArgs = formalArgs(IgEscapePlots), .cacheExtra = append(digestForPlot, digestOfData))
              },
              .args = argsAll)  # this tells mirai which objects are needed

            # IgEscapePlots(dt = data$covariates, bestModel = modelHere,
            #               climVar = sim$climateVariablesForFire[[igOrEsc]],
            #               rescalers = data$ignitionRescalers,
            #               fsProcess = "ignition", family = P(sim)$ignitionFamily,
            #               plotBiomass = P(sim)$plot_fuelBiomassPerPrediction,
            #               ignitionFitRTM = sim$ignitionFitRTM,
            #               studyAreaName = P(sim)$.studyAreaName,
            #               oPath = outputPath(sim)) |>
            #   Cache(omitArgs = c("dt", "bestModel"), .cacheExtra = digestOfData)
          }
          modelList
        })

  # put to sim sim$fireSense_IgnitionFitted, sim$fireSense_EscapeFitted
  list2env(fireSense_FittedModels, envir = envir(sim))
  on.exit() # let the mirai keep going, so remove the on.exit close of the mirai object
  return(invisible(sim))


  # ##############################
  # if ("ignition" %in% P(sim)$whichProcessesToFit) {
  #   digestIgnitionData <- .robustDigest(mget(c("fireSense_ignitionCovariates", "fireSense_ignitionFormula"),
  #                                            envir = envir(sim)))
  #
  #   ignitionData <- prepareCovariates(formula = sim$fireSense_ignitionFormula,
  #                                     covariates = sim$fireSense_ignitionCovariates,
  #                                     rescaleVars =P(sim)$rescaleVars) |>
  #     Cache(omitArgs = c("covariates", "formula"), .cacheExtra = digestIgnitionData)
  #
  #
  #   ignitionModel <- buildModel(covariates = ignitionData$covariates,
  #                               climVar = sim$climateVariablesForFire$ignition,
  #                               formula= ignitionData$formula, type = "ignition",
  #                               family = P(sim)$ignitionFamily) |>
  #     Cache(omitArgs = c("covariates", "formula"), .cacheExtra = digestIgnitionData)
  #
  #   #ignition specific
  #   origNoPix <- attributes(sim$ignitionFitRTM)$nonNAs   ## nrow(preSampleData) in eg above
  #   finalNoPix <- nrow(ignitionData$covariates)     ## nrow(postSampleData) in eg above
  #   lambdaRescaleFactor <- finalNoPix/origNoPix
  #
  #   modelList <- list(
  #     model = ignitionModel,
  #     rescales = ignitionData$ignitionRescalers,
  #     fittingRes = res(sim$ignitionFitRTM)[1],
  #     lambdaRescaleFactor = lambdaRescaleFactor)
  #
  #   sim$fireSense_IgnitionFitted <- modelList
  #   class(sim$fireSense_IgnitionFitted) <- "fireSense_IgnitionFit"
  #
  #   if (anyPlotting(P(sim)$.plots)) {
  #
  #     IgEscapePlots(dt = ignitionData$covariates, bestModel = ignitionModel,
  #                   climVar = sim$climateVariablesForFire$ignition,
  #                   rescalers = ignitionData$ignitionRescalers,
  #                   fsProcess = "ignition", family = P(sim)$ignitionFamily,
  #                   plotBiomass = P(sim)$plot_fuelBiomassPerPrediction,
  #                   ignitionFitRTM = sim$ignitionFitRTM,
  #                   studyAreaName = P(sim)$.studyAreaName,
  #                   oPath = outputPath(sim)) |>
  #       Cache(omitArgs = c("dt", "bestModel"), .cacheExtra = digestIgnitionData)
  #   }
  # }
  #
  # if ("escape" %in% P(sim)$whichProcessesToFit) {
  #   digestEscapeData <- .robustDigest(mget(c("fireSense_ignitionCovariates", "fireSense_ignitionFormula"),
  #                                            envir = envir(sim)))
  #
  #   escapeData <- prepareCovariates(formula = sim$fireSense_escapeFormula,
  #                                   covariates = sim$fireSense_escapeCovariates,
  #                                   rescaleVars = P(sim)$rescaleVars) |>
  #     Cache()
  #
  #   escapeModel <- buildModel(covariates = escapeData$covariates,
  #                             climVar = sim$climateVariablesForFire$ignition,
  #                             formula= escapeData$formula, type = "escape",
  #                             family = P(sim)$escapeFamily) |>
  #     Cache(omitArgs = c("covariates", "formula"), .cacheExtra = digestEscapeData)
  #
  #   modelList <- list(
  #     model = escapeModel,
  #     rescales = escapeData$ignitionRescalers,
  #     fittingRes = res(sim$ignitionFitRTM)[1],
  #     lambdaRescaleFactor = lambdaRescaleFactor)
  #
  #   sim$fireSense_EscapeFitted <- modelList
  #   class(sim$fireSense_EscapeFitted) <- "fireSense_EscapeFit"
  #
  #   if (anyPlotting(P(sim)$.plots)) {
  #     IgEscapePlots(dt = escapeData$covariates, bestModel = escapeModel,
  #                   climVar = sim$climateVariablesForFire$ignition,
  #                   fsProcess = "escape", family = P(sim)$escapeFamily,
  #                   rescalers = escapeData$ignitionRescalers,
  #                   plotBiomass = P(sim)$plot_fuelBiomassPerPrediction,
  #                   ignitionFitRTM = sim$ignitionFitRTM,
  #                   studyAreaName = P(sim)$.studyAreaName,
  #                   oPath = outputPath(sim)) |>
  #       Cache()
  #   }
  # }

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
      InterceptOnly <- grep(":", allTerms, value = TRUE, invert = TRUE)
      paste0(InterceptOnly, collapse = " + ")
    })
    paste0(allTermsNoMinus, collapse = " - ")
  })

  forms[["InterceptOnly"]] <- as.formula(paste0(terms[c(2,1,3)], collapse = " "), env = .GlobalEnv)
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

  dat <- m
  envir <- environment()
  system.time({
    mods <- Map(
      nam = names(forms), form = forms,
      # putting dat = m here causes it to become unresonsive --> a feature of "MoreArgs" in Map --> a list is evaluated
      MoreArgs = list(# dat = m, family = family,  # putting family = family here causes it to evaluated
                      type = type, envir = envir),
      f = function(form, nam, type, envir) {
        en <- new.env(parent = .GlobalEnv)
        if (type == "ignition") {
          ziform <- as.formula(paste0("~", paste0(climVar, collapse = "+")), env = en)
        } else {
          ziform <- as.formula(~0, env = en)
        }
        objNamesOutside <- c("dat", "family")
        objNamesInside <- c("form", "ziform", "nam")
        objsOutside <- mget(objNamesOutside, envir = envir)
        objsInside <- mget(objNamesInside)
        objs <- append(objsInside, objsOutside)
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
  whBest <- which.min(c(AICs[["full"]] + 2, AICs[["InterceptOnly"]], AICs[["climateOnly"]]))
  # whBest <- 1

  bestModel <- mods[[whBest]]
  # Remove the huge datasets
  bestModel$y <- NULL
  bestModel$data <- new.env(parent = emptyenv())
  messageColoured("Best model is:\n", messageFormulaFn(bestModel$call$formula), colour = "magenta")
  # summ <- summary(bestModel)

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


igOrEscNames <- function(igOrEsc, pre = "fireSense_", post, case = c("lower", "camel", "sentence", "title")) {
  if (startsWith(tolower(case[1]), prefix = "cam"))
    igOrEsc <- camelCase(igOrEsc)
  if (startsWith(tolower(case[1]), prefix = "sen") || startsWith(tolower(case[1]), prefix = "tit"))
    igOrEsc <- tools::toTitleCase(igOrEsc) # only has one word, so OK
  paste0(pre, igOrEsc, post)
}


camelCase <- function(x) {
  gsub("(^|[^[:alnum:]])([[:alnum:]])", "\\U\\2", x, perl = TRUE)
}
