defineModule(sim, list(
  name = "fireSense_IgnitionFit",
  description = paste("Fit statistical models that can be used to parameterize (calibrate)",
                      "the fire ignition component of landscape fire models (e.g. fireSense)."),
  keywords = c("fire frequency", "optimization", "additive property", "poisson",
               "negative binomial", "fireSense"),
  authors = c(
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut")),
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
                  "PredictiveEcology/fireSenseUtils@lccFix (>= 0.0.5.9065)",
                  "glmmTMB",
                  "ggplot2", "ggpubr", "MASS", "magrittr",
                  "numDeriv", "parallel", "parallelly",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/reproducible@AI (>= 2.1.2)",
                  #TODO correct this when reproducible is merged - it is due to cache(predict)
                  "RhpcBLASctl",
                  "PredictiveEcology/SpaDES.core@development (>= 2.0.2.9006)", "terra"),
  parameters = bindrows(
    defineParameter("family", "function, character", default = quote(poisson(link = "logit")),
                    desc = paste("a family function (must be wrapped with `quote()`) or a",
                                 "character string naming a family function.",
                                 "Only the negative binomial has been implemented",
                                 "For additional details see `?family`. This was formerly ",
                                 "`quote(MASS::negative.binomial(theta = 1, link = 'identity'))`.")),
    defineParameter("plot_fuelBiomassPerPrediction", "numeric", NULL, 1, 10,
                    desc = paste("when generating plots of climate x fuel class, the log of biomass (g/m2) for which",
                                 "to generate predictions across a gradient of climate values.",
                                 "If supplied, it will override any values in `sim$ignitionFitRTM$meanForestB`")),
    defineParameter("rescalers", "numeric", c("MDC" = 1000),
                    desc = paste("`NA` or a named vector of rescaling factors (numeric/integer values)",
                                 "for each predictor variable. If not `NA`, it will be used to rescale the",
                                 "variables as `var / rescalers['var']`. If `NA` and `rescaleVars == TRUE`,",
                                 "variables will be scaled to `[0,1]`.")),
    defineParameter("rescaleVars", "logical", default = FALSE,
                    desc = paste("Attempt to rescale variables? If `rescalers` is defined,",
                                 "use it to rescale variables as `var / rescalers['var']`. ",
                                 "Otherwise, `scale()` will be used to rescale variables to `[0,1]`,",
                                 "if they are not already within this range.")),
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
                 desc = paste("The column name(s) in the `fireSense_ignitionCovariates that is climate,",
                 "in a named list, .e.g. `climateVariablesForFire = list('ignition' = 'MDC')`")),
    expectsInput("fireSense_ignitionCovariates", "data.frame",
                 desc = "table of aggregated ignition covariates with annual ignitions"),
    expectsInput("flammableRTM", "SpatRaster", sourceURL = NA,
                 "RTM without ice/rocks/urban/water. Flammable map with 0 and 1."),
    expectsInput("ignitionFitRTM", "SpatRaster",
                 desc = paste("A (template) raster with information with regards to the spatial",
                              "resolution and geographical extent of `fireSense_ignitionCovariates.`",
                              "Used to pass this information onto `fireSense_ignitionFitted`",
                              "Needs to have number of non-NA cells as attribute:",
                              "(`ignitionFitRTM@data@attributes$nonNAs`), and optionally,",
                              "ignitionFitRTM@data@attributes$meanForestB")),
    expectsInput("fireSense_ignitionFormula", "character",
                 desc = "formula - as a character - describing the model to be fitted."),
  ),
  outputObjects = bindrows(
    createsOutput("covMinMax_ignition", "data.table",
                  desc = "Table of the original ranges (min and max) of covariates"),
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
  #TODO: review why Init event is not Init.
  if (!"pixelID" %in% colnames(sim$fireSense_ignitionCovariates)) {
    stop("fireSense_ignitionCovariates must have a 'pixelID' column")
  }

  if (is.empty.model(as.formula(sim$fireSense_ignitionFormula, env = .GlobalEnv))) {
    stop(moduleName, "> The formula describes an empty model.")
  }

  if (!all(is.na(P(sim)$rescalers))) {
    ## checks
    if (is.null(names(P(sim)$rescalers))) {
      stop("P(sim)$rescalers must be a named vector or NA.")
    }

    if (!all(names(P(sim)$rescalers) %in% names(sim$fireSense_ignitionCovariates))) {
      stop("names(P(sim)$rescalers) doesn't match variable names in fireSense_ignitionCovariates")
    }
  }

  if (is.null(attributes(sim$ignitionFitRTM)$nonNAs) ||
      length(attributes(sim$ignitionFitRTM)$nonNAs) == 0) {
    stop("sim$ignitionFitRTM@data@attributes$nonNAs must be a non-empty/non-NULL numeric")
  }

  return(invisible(sim))
}

frequencyFitRun <- function(sim) {

  moduleName <- current(sim)$moduleName

  fireSense_ignitionFormula <- as.formula(sim$fireSense_ignitionFormula, env = .GlobalEnv)
  terms <- terms.formula(fireSense_ignitionFormula)

  # Check the presence of at least one piecewise term
  fireSense_ignitionCovariates <- sim$fireSense_ignitionCovariates
  fireSense_ignitionCovariates <- copy(setDT(fireSense_ignitionCovariates))

  if (attr(terms, "response")) {
    y <- fireSense_ignitionFormula[[2L]]
  } else {
    stop(moduleName, "> Incomplete formula, the LHS is missing.")
  }

  ## any years in data?
  if (any(c("year", "yr") %in% tolower(names(fireSense_ignitionCovariates)))) {
    xvar <- intersect(c("year", "yr"), tolower(names(fireSense_ignitionCovariates)))
  } else {
    xvar <- rows
  }

  ## rescale variable and knots.
  if (isTRUE(P(sim)$rescaleVars)) {
    if (is.na(P(sim)$rescalers)) {
      ## TODO: lapply through each element in rescalers and assess which elements are to be rescaled vs normalized
      message("Variables outside of [0,1] range will be rescaled to [0,1]")

      needRescale <- fireSense_ignitionCovariates[, vapply(.SD, FUN = function(x) all(inrange(na.omit(x), 0, 1)),
                                                           FUN.VALUE = logical(1)),
                                                  .SDcols = notSpecialVars]
      message(paste("rescaling", needRescale))
      cols <- names(needRescale)[which(!needRescale)]
    } else {
      cols <- names(P(sim)$rescalers)
    }
    fireSense_ignitionCovariates <- rescaleVars(fireSense_ignitionCovariates, Par$rescalers)
  }

  # convert to data.table --> easier to work with
  m <- as.data.table(fireSense_ignitionCovariates)
  family <- P(sim)$family

  # Run 2 models ... the full one passed, plus a simplified one without interactions
  forms <- list()
  forms[["full"]] <- as.formula(fireSense_ignitionFormula, env = .GlobalEnv)
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

  # make minor modifications to dataset --> remove cases of >1 ignition per pixel (logit link for poisson model)
  climVar <- sim$climateVariablesForFire$ignition

  for (i in c("year", "ignitionsNoGT1"))
    set(m, NULL, i, as.integer(m[[i]]))
  for (i in c("yearChar"))
    set(m, NULL, i, factor(m[[i]]))

  # m <- m[sample(NROW(m), NROW(m)/4), ]
  system.time({
    mods <- Map(
      nam = names(forms), form = forms,
      MoreArgs = list(dat = m, family = family),
      function(form, nam, dat, family) {
        en <- new.env(parent = .GlobalEnv)
        ziform <- as.formula(paste0("~", paste0(climVar, collapse = "+")), env = en)
        # form <- as.formula(paste0("ignitionsNoGT1 ~ (1 | yearChar) + MDCc + youngAge + nonForest_highFlam + ",
        #                           "nonForest_lowFlam + class2 + class3"), env = en)
        objNames <- c("dat", "family", "form", "ziform", "nam")
        objs <- mget(objNames)
        dig <- en$dig <- .robustDigest(objs)
        list2env(objs, envir = en)
        message("Running glmmTMB with Zero-Inflated, Mixed effect, Poisson, using:\n",
                messageFormulaFn(form))

        out <- local({
          glmmTMB(form, data = dat,
                  ziformula = ziform, ## TODO this needs to be
                  family = eval(family)) |>
            ## Use .cacheExtra: there are lots of arguments to glmmTMB that seemed to be "always different"
            Cache(.functionName = paste0("glmmTMB_forIgnitions_", nam),
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
  whBest <- which.min(c(AICs[["full"]] + 2, AICs[["NoInteractions"]]))
  # whBest <- 1
  bestModel <- mods[[whBest]]
  messageColoured("Best model is:\n", messageFormulaFn(bestModel$call$formula), colour = "magenta")
  summ <- summary(bestModel)

  if (anyPlotting(P(sim)$.plots)) {
    ff <- as.character(bestModel$call$formula)

    formForNull <- as.formula(paste0(ff[[2]], ff[[1]], "1"), env = .GlobalEnv)

    nullModel <- glmmTMB(formForNull, dat = m, family = eval(family))

    ## https://stackoverflow.com/a/68684973 -- NOT VERY APPROPRIATE FOR RE model
    pseudoR2 <- as.numeric(1 - logLik(bestModel) / logLik(nullModel))

    message("Plotting has not been tested thoroughly")

    #### Plotting glmmTMB ####
    ## build two prediction datasets
    ## both predict across quantiles of climate variable
    ## the first dataset will have class means for cover and forest biomass
    ## the second will have cover values of 1 (i.e. representing complete cover for non-forest)
    ## and biomass values representing the mean for forested pixels (i.e. complete cover of forest)
    ## the 2nd dataset is different because the first one, while more representative of the actual landscape,
      ## implicitly includes non-forest pixels in the calculation of the mean.

    N <- 30
    p <- as.data.table(lapply(m, function(pp) if (is.numeric(pp)) mean(pp) else pp[1:N]))
    terms <- terms(bestModel)
    termsNonClimate <- setdiff(attr(terms, "term.labels"), climVar)

    ## populate a prediction dataset with quantiles of climate variable and mutually exclusive veg.
    for (var in climVar) {
      interpolateClimVar <- seq(quantile(m[[var]], 0.1),
                                (quantile(m[[var]], 0.95) * 1.5), length.out = N)
      set(p, NULL, var, interpolateClimVar)

      ## if more than 1 climate variable are used, they are plotted sequentially
      ## each prediction dataset will contain quantiles of one variable and mean of the other(s)
      otherClimVar <- climVar[!climVar %in% var]
      if (length(otherClimVar) > 0) {
        for (otherVar in otherClimVar) {
          set(p, NULL, otherVar, mean(m[[otherVar]]))
        }
      }

      termsNoInteraction <- termsNonClimate[termsNonClimate %in% names(m)]
      if (length(termsNoInteraction) == 0) { #all terms are interactions between fuel and climate
        termsWithInteraction <- attr(terms, "term.labels")
        termsNoInteraction <- sub(termsWithInteraction, pattern = climVar, replacement = "") |>
          sub(pattern = ":", replacement = "") #to catch ':<climvar>' or '<climVar>:'
      }

      #remove columns that aren't model terms - don't use set diff b/c some terms are only interaction
      unneededCovs <- setdiff(colnames(p), c(termsNoInteraction, climVar))
      for (rmCol in unneededCovs) {
        set(p, NULL, rmCol, NULL)
      }

      pAll <- rbindlist(lapply(seq(termsNoInteraction), function(x) p))
      pAll[, val := rep(termsNoInteraction, each = N)]

      termsUsingCover <- as.vector(m[, lapply(.SD, max), .SDcol = termsNoInteraction])
      termsUsingBiomass <- names(termsUsingCover[termsUsingCover > 1])
      termsUsingCover <- setdiff(names(termsUsingCover), termsUsingBiomass)

      for (val1 in c(termsUsingBiomass, termsUsingCover)) {
        set(pAll, which(!pAll$val %in% val1), val1, 0)
      }

      system.time({
        preds <- predict(object = bestModel, newdata = pAll, se.fit = TRUE, re.form = NA) |>
          Cache(omitArgs = "object", .cacheExtra = forms[whBest])
      })
      pAll[, pred := expit(preds$fit)]
      pAll[, val1 := factor(val)]
      pAll[, upper := expit(preds$fit + preds$se.fit)]
      pAll[, lower := expit(preds$fit - preds$se.fit)]

      resInKm2 <- prod(res(sim$ignitionFitRTM)) / 1e6 ## 1e6 m^2 == 1 km^2
      labelToUse <- paste("Ignition rate per", resInKm2, "km^2")
      filenameToUse <- paste0("IgnitionRatePer", resInKm2, "km2_", P(sim)$.studyAreaName, "_meanByClass")

      titl <- paste0("fireSense_IgnitionFit:", P(sim)$.studyAreaName,
                     " (", basename(outputPath(sim)), ")",
                     " -- Pseudo ")
      titl2 <- paste0(round(pseudoR2, 3))
      if (isTRUE(Par$rescaleVars)) {
        pAll <- rescaleVars(pAll, 1/Par$rescalers) # invert it
      }

      Plots(data = pAll, fn = plotFnLogitIgnition, # xColName = colName,
            ggylab = labelToUse,
            subtitle = "using mean cover and biomass per pixel",
            fillTitle = "veg. covariate",
            .plotInitialTime = NULL, # this means "ignore what `.plotInitialTime says; use only .plots`
            # centred = centred,
            climateVar = var, #TODO: fix to allow multiple climate variables
            # origXmax = max(sim$fireSense_ignitionCovariates[[colName]]), ## if supplied, adds bar to plot
            ggTitle = bquote(.(titl)~R^2 == .(titl2)),
            rawClimate =  m[[var]],
            filename = filenameToUse)

      ## make second prediction using mean forest or alternatively 100% non-forest cover

      if (!is.null(attributes(sim$ignitionFitRTM)$meanForestB) ||
          !is.null(P(sim)$plot_fuelBiomassPerPrediction)) {

        Bunit <- ifelse(!is.null(P(sim)$plot_fuelBiomassPerPrediction),
                        P(sim)$plot_fuelBiomassPerPrediction, #should be log already
                        log(attributes(sim$ignitionFitRTM)$meanForestB))
        BunitForLabel <- round(exp(Bunit), digits = 0)
        pAll2 <- copy(pAll)

        for (val1 in termsUsingBiomass) {
          set(pAll2, which(pAll2$val %in% val1), val1, Bunit)
        }

        for (val1 in termsUsingCover) {
          set(pAll2, which(pAll2$val %in% val1), val1, 1) #set the variable to 1 representing complete cover
        }

        system.time({
          preds <- predict(object= bestModel, newdata = pAll2, se.fit = TRUE, re.form = NA) |>
            Cache(omitArgs = "object", .cacheExtra = forms[whBest])
        })

        pAll2[, pred := expit(preds$fit)]
        pAll2[, val1 := factor(val)]
        pAll2[, upper := expit(preds$fit + preds$se.fit)]
        pAll2[, lower := expit(preds$fit - preds$se.fit)]

        filenameToUse <- paste0("IgnitionRatePer", resInKm2, "km2_", P(sim)$.studyAreaName, "_fullCoverAndBiomass")
        Plots(data = pAll2, fn = plotFnLogitIgnition,
              ggylab = labelToUse,
              subtitle = paste0("per ", BunitForLabel, " g B/m2 or 100% cover"),
              fillTitle = "veg. covariate",
              .plotInitialTime = NULL, # this means "ignore what `.plotInitialTime says; use only .plots`
              climateVar = var,
              rawClimate = m[[var]],
              # origXmax = max(sim$fireSense_ignitionCovariates[[colName]]), ## if supplied, adds bar to plot
              ggTitle = bquote(.(titl)~R^2 == .(titl2)),
              filename = filenameToUse)
      }
    }
    system.time({
      fittedNoRE <- predict(object = bestModel, newdata = m, se.fit = FALSE, re.form = NA,
                            type = "response") |>
        Cache(.functionName = "predict_forFitted_v_Obs_Ignitions",
              omitArgs = "object", .cacheExtra = forms[whBest])
    })

    # fittedVals <- fitted(bestModel)

    plotData <- data.table(fireSense_ignitionCovariates)
    plotData[,  rows := 1:nrow(plotData)]
    cols <- unique(c(paste(y), xvar, "rows"))
    plotData <- plotData[, ..cols]
    plotData <- cbind(plotData, fittedNoRE = fittedNoRE)

    predDT <- rbindlist(lapply(1:100,  FUN = function(x, DT) {
      rpoisPred <- rpois(nrow(DT), lambda = DT$fittedNoRE)
      n <- rep(x, nrow(DT))
      data.table(rpoisPred = rpoisPred, n = n, rows = DT$rows)
    }, DT = plotData))

    plotData <- plotData[predDT, on = "rows"]

    plotData <- plotData[, list(obsFires = sum(eval(y), na.rm = TRUE),
                                predFires = sum(rpoisPred, na.rm = TRUE)),
                         by = c(xvar, "n")]
    plotData[, obsFires := as.integer(obsFires)]
    plotData[, predFires := as.integer(predFires)]

    pd <- plotData[, .(obsFires = mean(obsFires), predFires = mean(predFires)), .(year)]
    correl <- cor(pd$obsFires, pd$predFires)

    plotData <- melt(plotData, id.var = c(xvar, "n"))

    Plots(data = plotData, fn = fittedVsObservedPlot,
          xColName = xvar, .plotInitialTime = NULL,
          ggylab = "num. fires",
          ggTitle = paste(P(sim)$.studyAreaName, "fireSense_IgnitionFit: obs. vs. fit"),
          ggSubtitle = paste0("Correlation = ", round(correl, 2)),
          filename = paste0("ignition_NumFiresFitted_", P(sim)$.studyAreaName))
  }

  mod$rescales <- if (isTRUE(P(sim)$rescaleVars)) {
    if (!all(is.na(P(sim)$rescalers))) {
      ## TODO: allow list of rescalers to be passed with mix of NA and other vals
      sapply(needRescale, FUN = function(x) {
        paste0("fireSenseUtils::rescale(", x, ", to = c(0,1))")
      }, USE.NAMES = TRUE, simplify = FALSE)
    } else {
      sapply(names(P(sim)$rescalers), FUN = function(x, vec) {
        paste(x, "/", vec[x])
      }, vec = P(sim)$rescalers, USE.NAMES = TRUE, simplify = FALSE)
    }
  } else {
    NULL
  }

  origNoPix <- attributes(sim$ignitionFitRTM)$nonNAs   ## nrow(preSampleData) in eg above
  finalNoPix <- nrow(fireSense_ignitionCovariates)     ## nrow(postSampleData) in eg above
  lambdaRescaleFactor <- finalNoPix/origNoPix

  modelList <- list(
    model = bestModel,
    #formula, data, coef, coef.se, convergence should all be attainable.
    # formula = forms[[whBest]],
    # convergence = bestModel$fit$convergence,
    rescales = mod$rescales,
    fittingRes = res(sim$ignitionFitRTM)[1],
    lambdaRescaleFactor = lambdaRescaleFactor)

  sim$fireSense_IgnitionFitted <- modelList
  class(sim$fireSense_IgnitionFitted) <- "fireSense_IgnitionFit"

  return(invisible(sim))
}

fittedVsObservedPlot <- function(d, ggTitle, ggSubtitle = NULL, ggylab, xColName)  {
  ggplot <- ggplot(data = d, aes_string(x = xColName, y = "value", colour = "variable")) +
    stat_summary(aes(fill = variable), fun.data = mean_ci,
                 geom = "ribbon", alpha = 0.5, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    scale_color_discrete(labels = c("obsFires" = "observed no. fires",
                                    "predFires" = "fitted no. fires")) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(y = ggylab, x = xColName, title = ggTitle,
         subtitle = ggSubtitle, colour = "")
  ggplot
}

expit <- function(x) 1/(1 + exp(-x)) # inverse logit function; used below

messageFormulaFn <- function(form) {
  gsub(" {2,100}", " ", paste0(format(form), collapse = ""))
}

plotFnLogitIgnition <- function(pAll, subtitle = NULL, ggylab,
                                ggTitle, fillTitle, climateVar, rawClimate) {

  quants <- quantile(rawClimate, probs = c(0.25, 0.5, 0.75, 0.95))
  ggplot(pAll, aes(x = .data[[climateVar]], # MDCc + centred,
                   y = pred, by = val1, col = val1)) +
    geom_line(aes(y = pred), lwd = 1.5) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = val1), alpha = 0.4, show.legend = FALSE) +
    geom_vline(xintercept = quants, show.legend = FALSE, linetype = "dotted") +
    annotate("text", y = max(pAll$pred),
             x = quants - mean(quants*0.025),
             label = c("25%", "50%", "75%", "95%"), angle = 90) +
    labs(y = ggylab, title = ggTitle, col = fillTitle, subtitle = subtitle) +
    guides(lwd = "none", alpha = "none") +
    theme_bw()
}

#TODO: this function is currently unused - remove it if no longer needed
# identifyEnvs <- function(l, topEnv) {
#   if (is.list(l)) {
#     out <- lapply(l, function(ll) {
#       identifyEnvs(ll, topEnv)
#     })
#   } else {
#     if (NROW(l) == 0 || is(l, "externalptr")) {
#       return(NULL)
#     } else if (is.function(l)) {
#       return(environment(l))
#     } else if (is.environment(l)) {
#       return(l)
#     }
#     return(NULL)
#   }
#   return(out)
# }

rescaleVars <- function(dt, rescalers) {
  cols <- names(Par$rescalers)
  dt[, (cols) := mapply(FUN = function(x, vec) {x / vec},
                        x = .SD, vec = rescalers,
                        SIMPLIFY = FALSE),
     .SDcols = cols]
  dt[]
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
