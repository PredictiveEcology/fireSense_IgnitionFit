defineModule(sim, list(
  name = "fireSense_IgnitionFit",
  description = paste("Fit statistical models that can be used to parameterize (calibrate)",
                      "the fire ignition component of landscape fire models (e.g. fireSense)."),
  keywords = c("fire frequency", "optimization", "additive property", "poisson",
               "negative binomial", "fireSense"),
  authors = c(
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "aut"),
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = "aut"),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = character(),
  version = list(fireSense_IgnitionFit = "1.0.1.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_IgnitionFit.Rmd"),
  loadOrder = list(after = "fireSense_dataPrepFit",
                   before = "fireSense_dataPrepPredict"),
  reqdPkgs = list("data.table", "dplyr", "PredictiveEcology/SpaDES.core@development (>= 3.0.4)",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.1.0)",
                  "glmmTMB",
                  "ggplot2", "ggpubr", "MASS", "magrittr",
                  "numDeriv", "parallel", "parallelly",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/reproducible@development (>= 2.1.2.9067)",
                  #TODO correct this when reproducible is merged - it is due to cache(predict)
                  "RhpcBLASctl",
                  "caret", "pROC",
                  "PredictiveEcology/SHAPforxgboost (>= 0.1.3.9001)", "xgboost (>=3.0.0)", "lightgbm", # install.packages('xgboost', repos = c('https://dmlc.r-universe.dev', 'https://cloud.r-project.org'))
                  "terra"),
  parameters = bindrows(
    defineParameter("crossValType", "character", c("time-ordered", "crossValidation"), NA, NA,
                    paste("Currently has no effect on the fit: `runXGBOOST()` always uses k-fold cross validation",
                          "(its time-ordered branch is disabled). The first element is used in the plot filename.")),
    defineParameter("escapeFamily", "function, character", default = quote(binomial(link = "logit")),
                    desc = paste("Currently unused. A family function (wrapped with `quote()`) or a",
                                 "character string naming one, for the escape model of non-xgboost algorithms.")),
    defineParameter("ignitionFamily", "function, character", default = quote(poisson(link = "log")),
                    desc = paste("Currently unused. A family function (wrapped with `quote()`) or a",
                                 "character string naming one, for the ignition model of non-xgboost algorithms.")),
    defineParameter("plot_fuelBiomassPerPrediction", "numeric", NULL, 1, 10,
                    desc = "Currently unused."),
    defineParameter("rescaleVars", "logical", default = TRUE,
                    desc = paste("If `TRUE`, covariates are centred and scaled with `scale()` before fitting.",
                                 "The centring and scaling values are returned in the `scaleData` element of the outputs.")),
    defineParameter("whichProcessesToFit", "character", c("ignition", "escape"), NA, NA,
                    "which processes to fit: ignition, escape, or both (the default)"),
    defineParameter("modelAlgorithm", "character", "xgboost", NA, NA,
                    "Model type. Only `xgboost` (any value containing 'xgb') works."),
    defineParameter(".plots", "character", default = "screen",
                    desc = paste("See `?Plots`. If set, plots the predicted response against each covariate,",
                                 "for climate and fuel covariates separately, and saves it as png in `figurePath(sim)`.")),
    defineParameter(".plotInitialTime", "numeric", default = NULL,
                    desc = "Unused. Plots are made in the `run` event if `.plots` is set."),
    defineParameter(".runInitialTime", "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start time of the simulation."),
    defineParameter(".runInterval", "numeric", default = NA,
                    desc = paste("optional. Interval between two runs of this module,",
                                 "expressed in units of simulation time.",
                                 "By default, NA, which means that this module only runs once per simulation.")),
    defineParameter(".saveInitialTime", "numeric", default = NA,
                    desc = "Leave as `NA`: the module has no `save` event, so a value only produces a warning."),
    defineParameter(".saveInterval", "numeric", default = NA,
                    desc = "Unused."),
    defineParameter(".seed", "list", NULL, NA, NA,
                    paste("Named list of seeds to use for each event (names).",
                          "E.g., `list('init' = 123)` will `set.seed(123)`",
                          "at the start of the init event and unset it at the end.",
                          "Defaults to `NULL`, meaning that no seeds will be set.")),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Currently unused."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules,",
                                 "where stochasticity and time are not relevant."))
  ),
  inputObjects = bindrows(
    expectsInput("climateVariablesForFire", "list",
                 desc = paste("Named list of the climate column name(s) in `fireSense_ignitionCovariates`,",
                              "e.g. `list('ignition' = 'MDC')`. Not used by xgboost.")),
    expectsInput("fireSense_ignitionCovariates", "data.frame",
                 desc = paste("Table of aggregated ignition covariates with annual `ignitions` counts,",
                              "one row per `pixelID` and `year`.")),
    expectsInput("ignitionFitRTM", "SpatRaster",
                 desc = paste("Template raster at the resolution and extent of `fireSense_ignitionCovariates`.",
                              "Its resolution and its number of non-NA cells, which must be in the attribute `nonNAs`",
                              "(`attributes(ignitionFitRTM)$nonNAs`), are stored in the fitted objects.")),
    expectsInput("fireSense_ignitionFormula", "character",
                 desc = "Model formula, as character. Not used by xgboost."),
  ),
  outputObjects = bindrows(
    createsOutput("fireSense_EscapeFitted", "fireSense_EscapeFit",
                  desc = paste("List of `modelList` and `scaleData`, as `fireSense_IgnitionFitted`,",
                               "with `modelList` of class `fireSense_EscapeFit`.")),
    createsOutput("fireSense_IgnitionFitted", "fireSense_IgnitionFit",
                  desc = paste("List of `modelList` (class `fireSense_IgnitionFit`: `model`, the per-fold xgboost models",
                               "and their ROC curves, `fittingRes`, `lambdaRescaleFactor`, `rescales`, `family`)",
                               "and `scaleData` (centre and scale used to standardise the covariates)."))
  )
))

## event types
#   - type `init` is required for initialiazation

#' Event dispatcher
#'
#' `init` schedules `checkData` and `run` at `.runInitialTime`; `run` repeats every
#' `.runInterval` if that is not `NA`.
#'
#' @param sim A `simList`.
#' @param eventTime Time of the event.
#' @param eventType One of `init`, `checkData`, `run`.
#' @param debug Unused.
#' @return `sim`, invisibly.
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
#' Check inputs before fitting (the `checkData` event)
#'
#' Requires the `nonNAs` attribute on `ignitionFitRTM` when fitting ignition.
#' For non-xgboost algorithms also runs `checkData()`.
#'
#' @param sim A `simList`.
#' @return `sim`, invisibly, unchanged.
Init <- function(sim) {

  if ("ignition" %in% P(sim)$whichProcessesToFit) {
    if (is.null(attributes(sim$ignitionFitRTM)$nonNAs) ||
        length(attributes(sim$ignitionFitRTM)$nonNAs) == 0) {
      stop("sim$ignitionFitRTM@data@attributes$nonNAs must be a non-empty/non-NULL numeric")
    }
    if (grepl("xgb", Par$modelAlgorithm) %in% FALSE)
      checkData(formula = sim$fireSense_ignitionFormula,
                covariates = sim$fireSense_ignitionCovariates, type = "ignition")
  }
  if ("escape" %in% P(sim)$whichProcessesToFit) {
    if (grepl("xgb", Par$modelAlgorithm) %in% FALSE)
      checkData(sim$fireSense_escapeFormula, sim$fireSense_escapeCovariates, "escape")
  }
  if (!any(c("ignition", "escape") %in% P(sim)$whichProcessesToFit)) {
    stop("please review P(sim)$whichProcesesToFit...ensure lower-case")
  }

  return(invisible(sim))
}

#' Fit the requested processes (the `run` event)
#'
#' Calls `buildModelsFitModels()` for each of `P(sim)$whichProcessesToFit` and assigns the
#' results to `sim$fireSense_IgnitionFitted` and/or `sim$fireSense_EscapeFitted`.
#'
#' @param sim A `simList`.
#' @return `sim`, invisibly.
frequencyFitRun <- function(sim) {

  #TODO: make cache smart by digest args in advance
  whichProcessesToFit <- P(sim)$whichProcessesToFit
  names(whichProcessesToFit) <- igOrEscNames(whichProcessesToFit, post = "Fitted", case = "sen")

  fireSense_FittedModels <-
    purrr::pmap(.l  = list(igOrEsc = whichProcessesToFit), sim = sim,
                .f = buildModelsFitModels)

  # put to sim sim$fireSense_IgnitionFitted, sim$fireSense_EscapeFitted
  list2env(fireSense_FittedModels, envir = envir(sim))
  return(invisible(sim))
}

#' Formula as a one-line string
#'
#' @param form A formula.
#' @return Character string, with runs of spaces collapsed.
messageFormulaFn <- function(form) {
  gsub(" {2,100}", " ", paste0(format(form), collapse = ""))
}

#' Validate covariates and formula (non-xgboost algorithms only)
#'
#' @param formula Model formula, or character coercible to one.
#' @param covariates data.frame; must have a `pixelID` column.
#' @param type `"ignition"` or `"escape"`; used in the error messages.
#' @return `NULL`, invisibly. Stops if a check fails.
checkData <- function(formula, covariates, type) {
  if (!"pixelID" %in% colnames(covariates)) {
    stop("covariates for ", type, " must have a 'pixelID' column")
  }

  if (is.empty.model(as.formula(formula, env = .GlobalEnv))) {
    stop("formula for ", type, " describes an empty model.")
  }
}


#' Fit the model for one process
#'
#' With an xgboost `modelAlgorithm` this only calls `runXGBOOST()`. Other algorithms fit the
#' formulas in `formsToRun` and keep the one with lowest AIC (the full model penalised by 2),
#' but that path cannot currently be reached: `fireSenseUtils::prepareCovariatesOuter()`
#' stops for them.
#'
#' @param covariates data.table of response and covariates. The response columns are
#'   converted to integer by reference.
#' @param formula Model formula. Not used by xgboost.
#' @param type `"ignition"` or `"escape"`.
#' @param climVar Climate column name. Not used by xgboost.
#' @param family Model family. Not used by xgboost.
#' @param digestOfData List of digests; its `fireSense_ignitionCovariates` element is passed
#'   to `runXGBOOST()` as `dig`.
#' @param modelAlgorithm See the `modelAlgorithm` module parameter.
#' @param nFolds Number of cross-validation folds.
#' @param crossValType Passed to `runXGBOOST()`, which ignores it.
#' @param formsToRun Which formula variants to fit. Not used by xgboost.
#' @return For xgboost, see `runXGBOOST()`. Otherwise the best fitted model with its data removed.
buildModel <- function(covariates, formula, type = "ignition",
                       climVar,
                       family,
                       digestOfData,
                       modelAlgorithm, nFolds = 5, crossValType = c("time-ordered", "crossValidation"),
                       formsToRun = c("full", "InterceptOnly", "climateOnly")) {

  dat <- covariates

  # Run 3 models ... the full model, climate only,  and NULL model
  if (any(grepl("xgb", modelAlgorithm) %in% FALSE)) {
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
  }
  varToInt <- fireSenseUtils::ignitionsTxt
  if (type == "escape") {
    varToInt <- c(fireSenseUtils::ignitionsTxt, fireSenseUtils::escapesTxt)
  }
  if (!any(grepl("xgb", modelAlgorithm)))
    varToInt <- c(varToInt, "year")

  for (i in varToInt)
    set(dat, NULL, i, as.integer(dat[[i]]))

  envir <- environment()
  mods <- list()

  if (any(grepl("xgb", modelAlgorithm))) {
    nam <- type
    mods[[nam]] <- runXGBOOST(dat, digestOfData$fireSense_ignitionCovariates,
                              type = type, nFolds = nFolds, crossValType = crossValType[1])
  } else {
    forms <- forms[formsToRun]

    for (i in c("yearChar"))
      set(dat, NULL, i, factor(dat[[i]]))
    for (ind in seq_along(forms)) {
      nam <- names(forms)[[ind]]
      form <- forms[[ind]]
      if (any(grepl("adaptive", modelAlgorithm))) {
        mods[[nam]] <- runGLMMAdaptiveWithSimplifications(nam, form, dat, family, type, climVar)
      } else if (any(grepl("nb", modelAlgorithm))) {
        mods[[nam]] <- runGLM.NB(dat)
      } else {
        stop("Other modelAlgorithms not implemented")
      }
    }
  }

  # AIC is not relevant for the xgboost
  if (any(grepl("xgb", modelAlgorithm) %in% FALSE)) {
    AICs <- sapply(mods, AIC)
    ## even if the AIC is <2 better, should take simpler model;
    ## in tests, turned many to non-significant when had interactions
    whBest <- which.min(c(AICs[["full"]] + 2, AICs[["InterceptOnly"]], AICs[["climateOnly"]]))

    bestModel <- mods[[whBest]]
    # Remove the huge datasets
    bestModel$y <- NULL
    bestModel$data <- new.env(parent = emptyenv())
    messageColoured("Best model is:\n", messageFormulaFn(bestModel$call$formula), colour = "magenta")
  } else {
    bestModel <- mods[[type]]
  }

  return(bestModel)
}


#' Default inputs
#'
#' Stops if `fireSense_ignitionCovariates` is not supplied. If `climateVariablesForFire` is
#' not supplied, sets it to `list(ignition = "MDC")`, which requires an `MDC` covariate column.
#'
#' @param sim A `simList`.
#' @return `sim`.
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

#' Prepare covariates, fit one process, and plot it
#'
#' @param igOrEsc `"ignition"` or `"escape"`. Covariates are read from
#'   `sim$fireSense_ignitionCovariates` or `sim$fireSense_escapeCovariates`.
#' @param sim A `simList`.
#' @return List of `modelList`, of class `fireSense_IgnitionFit` or `fireSense_EscapeFit`, and
#'   `scaleData`, the centre and scale used to standardise the covariates. `modelList` holds
#'   `model` (from `buildModel()`), `rescales`, `fittingRes` (resolution of `ignitionFitRTM`),
#'   `lambdaRescaleFactor` (rows of covariates / non-NA cells of `ignitionFitRTM`) and `family`.
buildModelsFitModels <- function(igOrEsc, sim) {
  covariatesHere <- igOrEscNames(igOrEsc, post = "Covariates")
  data <- prepareCovariatesOuter(sim[[covariatesHere]], 
                                 algorithm = Par$modelAlgorithm, 
                                 rescaleVars = Par$rescaleVars)
  
  nFolds <- 5

  crossValType <- Par$crossValType

  modelHere <- buildModel(covariates = data$covariates,
                          climVar = sim$climateVariablesForFire$ignition,
                          formula= data$formula,
                          type = igOrEsc, nFolds = nFolds,
                          crossValType = crossValType,
                          family = family,
                          digestOfData = data$digestOfData,
                          modelAlgorithm = Par$modelAlgorithm,
                          formsToRun = c("full", "InterceptOnly", "climateOnly")[1]) |>
    Cache(omitArgs = c("covariates", "formula"),
          .functionName = paste0("buildModel for ", igOrEsc),
          .cacheExtra = c(data$digestOfData, 1), cacheSaveFormat = "rds") # doesn't work with qs
  #ignition specific
  origNoPix <- attributes(sim$ignitionFitRTM)$nonNAs   ## nrow(preSampleData) in eg above
  finalNoPix <- nrow(data$covariates)     ## nrow(postSampleData) in eg above
  lambdaRescaleFactor <- finalNoPix/origNoPix

  modelList <- list(
    model = modelHere,
    rescales = data$ignitionRescalers,
    fittingRes = res(sim$ignitionFitRTM)[1],
    lambdaRescaleFactor = lambdaRescaleFactor,
    family = family)

  class(modelList) <- igOrEscNames(igOrEsc, post = "Fit", case = "Title")

  if (anyPlotting(P(sim)$.plots)) {
    message("Plotting ", igOrEsc, "...")

    figPath <- figurePath(sim)
    digModels <- attr(modelHere, "tags")
    modelOnly <- modelHere[grep("Fold", names(modelHere))]
    aa <- setupPlots(modelOnly, dat = data$covariates, igOrEsc) |>
      Cache(omitArgs = c("dat", "modelOnly"),
            .cacheExtra = list(digestOfData = data$digestOfData, digModels = digModels, plotPredictions = plotPredictions),
            .functionName = paste0(".functionName_", igOrEsc))
    fn1 <- functionNameHelper("FuelClimate", ifelse(igOrEsc == "escape", "", "Lightning"), "predicted", igOrEsc, Par$crossValType[1])
    fn <- functionNameHelper(fn1, format(Sys.time()))
    bb <- quote(ggarrange(plotlist = aa))

    a <- Plots(bb, path = figPath,
               filename = fn,
               type = "png", ggsaveArgs = list(width = 8, height = 5, scale = 1.5), useCache = TRUE)
  }
  list(modelList = modelList, scaleData = attr(data$covariates, "scaleData"))
}


#' Fit a zero-inflated mixed model, simplifying it on failure
#'
#' Up to three attempts: `GLMMadaptive::mixed_model()` with zi.negative.binomial, then
#' zi.poisson, then `lme4::glmer()` with poisson; or zi.binomial, then `glmer()` with binomial.
#' Cannot currently be reached; see `buildModel()`.
#'
#' @param nam Label of the formula variant; used in the cache function name.
#' @param form Formula with one random-effect term.
#' @param dat data.table of response and covariates.
#' @param family Quoted family call, e.g. `quote(poisson(link = "log"))`.
#' @param type `"ignition"` (zero inflation on `climVar`) or `"escape"` (none).
#' @param climVar Climate column name(s).
#' @return List of the attempted fits; the last element is the accepted one, or a failure.
runGLMMAdaptiveWithSimplifications <- function(nam, form, dat, family, type, climVar) {
  en <- new.env(parent = asNamespace("fireSenseUtils"))
  enObjs <- new.env(parent = en)
  if (type == "ignition") {
    ziform <- paste0("~", paste0(climVar, collapse = "+"))
  } else {
    ziform <- "~0"
  }
  objsNames <- c("form", "ziform", "family", "nam", "ranForm2", "form2", "thisFamily")

  other <- c("dat", "dat")

  ziform <- as.formula(ziform, env = .GlobalEnv)
  tt <- terms(form)
  whRanTerms <- grep("\\|", attributes(tt)$term.labels, value = TRUE)
  whRanTerms2 <- gsub("\\|", "\\\\|", whRanTerms)
  fixedTerms <- gsub(paste0("\\( *", whRanTerms2, " *\\) *\\+*"), "", form[-(1:2)])
  if (nzchar(fixedTerms) %in% FALSE)
    fixedTerms <- 0
  form2 <- as.formula(paste0(format(form[[2]]), format(form[[1]]), fixedTerms), env = .GlobalEnv)
  ranForm2 <- as.formula(paste0("~", whRanTerms), env = .GlobalEnv)
  ziBinomial <- "zi.binomial"
  ziNegativeBinomial <- "zi.negative.binomial"
  ziPoisson <- "zi.poisson"
  binomi <- "binomial"
  poisso <- "poisson"
  familyExtra <- data.table(orig = c(poisso, binomi), forMixedModel = list(ziNegativeBinomial,
                                                                           ziBinomial))
  fam <- format(family[[1]])
  thisFamily <- familyExtra[ orig %in% fam]$forMixedModel[[1]]

  objs <- mget(objsNames)
  model <- quote(GLMMadaptive::mixed_model(form2, random = ranForm2, data = dat, family = thisFamily,
                                           zi_fixed = ziform,
                                           control = list(iter_EM = 0, max_coef_value = 1000)))
  objs[["model"]] <- model
  dig <- .robustDigest(objs)
  other <- mget(other)
  dig <- append(dig, .robustDigest(other))

  objs[["dig"]] <- dig # put in other environment

  list2env(other, envir = en)

  out <- list()
  for (i in 1:3) {
    # Put the objects and functions in the other environnemt
    list2env(objs, envir = enObjs)

    oo <- capture.output(objs$model)
    message("Running ", oo, "\n", "where 'form' is: ", messageFormulaFn(form))

    library(GLMMadaptive)
    print(system.time(
      out[[i]] <- try(local(
        eval(model) |>
          Cache(.functionName = paste0("mixed_model_for", "_", nam),
                omitArgs = c(formalArgs(GLMMadaptive::mixed_model), "level", "x", "origDat"),
                .cacheExtra = dig)
        , envir = enObjs)
      )))


    ss <- try(summary(out[[i]]))
    # Lots of errors:
    # 1. Error in solve.default(object$Hessian) :
    #    system is computationally singular: reciprocal condition number = 2.30595e-131 --> tested with !is(ss, "try-error")
    # 2. Hessian matrix at convergence is not positive definite; unstable solution. --> tested with !is.na(ss$phis_table[[2]])

    if (isTRUE(
      !is(out[[i]], "try-error") && !is(ss, "try-error") &&
      tryCatch(!anyNA(ss$coef_table), error = function(e) TRUE) &&
      anyNA(ss$phis_table[[2]]) %in% FALSE)
    ) {
      break
    }
    useGLMER <- FALSE
    # if fails --> try simpler models zi.negative.binomial --> zi.poisson --> poisson
    # if fails --> try simpler models zi.binomial --> binomial
    if (identical(i, 1L)) {
      if (identical(objs$thisFamily, ziNegativeBinomial)) {
        objs$thisFamily <- ziPoisson
      }

      if (identical(objs$thisFamily, ziBinomial)) {
        objs$thisFamily <- binomi
        useGLMER <- TRUE
      }
    } else if (identical(i, 2L)) {
      if (identical(objs$thisFamily, ziPoisson)) {
        objs$thisFamily <- poisso
        useGLMER <- TRUE
      }

    }
    if (useGLMER) {
      glmerModel <- quote(lme4::glmer(form, data = dat, family = thisFamily))
      objs$model <- model <- glmerModel
      objs$dig$model <- .robustDigest(objs$model)
    }

    objs$dig$thisFamily <- .robustDigest(objs$thisFamily)

  }
  out
}

#' Fit cross-validated xgboost models
#'
#' One model per fold, Tweedie objective, trained on all rows except the fold's, which are the
#' `eval_set`. Folds come from `cvFolds()`. For `type = "escape"` only rows with
#' `ignitions > 0` are used. The RNG seed is fixed, and restored on exit, so folds repeat and
#' caching works.
#'
#' @param dat data.table of response and covariates. `pixelID`, `year` and `yearChar` columns
#'   are dropped.
#' @param dig Digest of `dat`, used to cache the creation of dummy variables.
#' @param type `"ignition"` or `"escape"`; also the pattern that finds the response column.
#' @param nFolds Number of folds.
#' @param crossValType Ignored. The time-ordered branch is disabled, so folds are always k-fold.
#' @return List of the `nFolds` xgboost models, named `Fold1`, ..., plus `rocs` from `rocPerFold()`.
runXGBOOST <- function(dat, dig, type = "ignition", nFolds = 5,
                       crossValType = c("time-ordered", "crossValidation")) {
  dat3Forxgboost <- dat

  # Add dummy variables for factor columns -- i.e., the random effects
  if (all(sapply(dat3Forxgboost, is.numeric)) %in% FALSE)
    dat3Forxgboost <-
      model.matrix(~ . + 0, data = dat3Forxgboost) |>
      Cache(omitArgs = c("object", "data", "x"), .cacheExtra = dig)# Creates dummy variables


  savedSeed <- .Random.seed
  on.exit(assign(".Random.seed", savedSeed, envir = .GlobalEnv), add = TRUE)
  set.seed(12345) # so kfolds are same, so Caching works correctly below; if dat3 changes number of rows,
  # it will be a totally different sequence; but it will be the same sequence
  # if number of rows doesn't change

  # Setup k-folds
  yearColname <- grep("year", colnames(dat3Forxgboost), value = TRUE)
  indexNames <- c("keepAll", "keepEval")

  if (identical(type, "escape")) {
    dat3Forxgboost <- dat3Forxgboost[ignitions > 0]
  }
  if (length(yearColname) && FALSE) {
    crossValType <- "time-ordered"
    times <- unique(dat3Forxgboost[[yearColname]])
    testLength <- 3
    initialWindow <- length(times) - testLength - nFolds + 1
    trainIndexK <- caret::createTimeSlices(times, initialWindow = initialWindow, testLength, fixedWindow = FALSE)
    trainIndexK <- Map(tr = trainIndexK$train, te = trainIndexK$test, function(tr, te) {
      keepAll <- which(dat3Forxgboost[[yearColname]] %in% times[c(tr, te)])
      keepEval <- which(dat3Forxgboost[[yearColname]] %in% times[te])
      list(keepAll, keepEval) |> setNames(indexNames)
    })

  } else {
    crossValType <- "crossValidation"
    trainIndexK <- cvFolds(dat3Forxgboost[[grep(value = TRUE, type, colnames(dat3Forxgboost))]], nFolds)
  }

  colOrder <- setdiff(colnames(dat3Forxgboost), c("pixelID", "year"))
  colOrder <- colOrder[colOrder %in% grep("yearChar", colnames(dat3Forxgboost), invert = TRUE, value = TRUE)]
  colOrder <- sample(colOrder)
  dat3Forxgboost <- dat3Forxgboost[, ..colOrder]
  dig <- .robustDigest(dat3Forxgboost)
  colnamesNoIgn <- grep(paste0(fireSenseUtils::ignitionsTxt,"|",fireSenseUtils::escapesTxt), colnames(dat3Forxgboost), value = TRUE, invert = TRUE) |>
    sort() # make alphabetical
  dat3ForxgboostNoIgn <- dat3Forxgboost[, ..colnamesNoIgn]

  ignOrEscapeColName <- grep(value = TRUE,type, colnames(dat3Forxgboost))
  stopIfFoldsLackPositives(dat3Forxgboost[[ignOrEscapeColName]], trainIndexK, type)
  st <- system.time(
    mm <- purrr::pmap(
      list(valInd = trainIndexK, kFold = seq(nFolds)),
      function(valInd, kFold)
      {

        wholeDataset <- valInd[[indexNames[[1]]]]
        valInd <- valInd[[indexNames[[2]]]]
        digValInd <- .robustDigest(valInd) # should be eval set

        # xgboost objects do not save with `qs` ... must be `rds`
        mTweedie <- xgboost(x = dat3ForxgboostNoIgn[wholeDataset], # should be whole set
                            y = dat3Forxgboost[, get(ignOrEscapeColName)][wholeDataset],
                            objective = "reg:tweedie",
                            nthread = 10,
                            eval_set = valInd, # should be keepEval
                            monitor_training = TRUE,
                            verbosity = 0,
                            nrounds = 100,
                            max_depth = 2,
                            reg_lambda = 0.5,
                            learning_rate = 0.15
        ) |> Cache(omitArgs = c("x", "y", "eval_set"),
                   .functionName = functionNameHelper("xgboost", type, kFold),
                   .cacheExtra = c(dig, digValInd, type),
                   cacheSaveFormat = "rds")
        # Predict probabilities
        valData <- dat3Forxgboost[valInd, ]
        pred2 <- predict(mTweedie, valData)
        if (!all(is.finite(pred2))) {
          stop(type, " cross-validation fold ", kFold, ": the model predicted ", sum(!is.finite(pred2)),
               " non-finite values out of ", length(pred2), call. = FALSE)
        }
        valData <- cbind(valData, predTweedie = pred2)

        shap_values <- shap.values(mTweedie, dat3ForxgboostNoIgn) |>
          Cache(omitArgs = formalArgs(shap.values),
                .functionName = functionNameHelper("shap.values", type, kFold),
                .cacheExtra = c(dig, digValInd, type))
        shapContrib <- shap_values$shap_score
        shapContrib <- shapContrib[, -"(Intercept)"]
        shap_long <- shap.prep(
          shap_contrib = shapContrib, X_train = dat3ForxgboostNoIgn)  |>
          Cache(omitArgs = formalArgs(shap.prep),
                .functionName = functionNameHelper("shap.prep", type, kFold),
                .cacheExtra = c(dig, digValInd, type))
        list(valData = valData, mod = mTweedie, shap_long = shap_long)
      })
  )
  rocs <- rocPerFold(mm, ignOrEscapeColName)
  tweedie <- aucPerFold(rocs)
  print(meanRocMessage(tweedie))
  mm2 <- Map(m = mm, function(m) m$mod)
  mm2 <- append(mm2, list(rocs = rocs))

  return(mm2)
}


#' ROC curve per validation fold
#'
#' AUC needs both outcomes present. A validation fold can legitimately hold only one --
#' for escape in a small study area, every fire escaped or none did -- and `pROC::roc()`
#' stops there with "'response' must have two levels". That used to abort the whole fit,
#' which is the wrong trade: these curves are a diagnostic, printed and returned alongside
#' the models, and nothing about the fitted model (`m$mod`) depends on them. A fold with
#' one outcome therefore reports no AUC and the fit stands.
#'
#' @param mm List with one element per fold, each holding `valData`: the validation rows with
#'   the response and the prediction, `predTweedie`.
#' @param ignOrEscapeColName Name of the response column.
#' @return List of `pROC::roc` objects; `NULL` for a fold with only one outcome.
rocPerFold <- function(mm, ignOrEscapeColName) {
  lapply(mm, function(d) {
    ignZeroAndOnes <- pmin(1L, d$valData[[ignOrEscapeColName]])
    if (length(unique(ignZeroAndOnes)) < 2L) {
      message("  ... a validation fold of ", ignOrEscapeColName, " holds only the value ",
              unique(ignZeroAndOnes), "; AUC is undefined for it, so it is skipped")
      return(NULL)
    }
    pROC::roc(ignZeroAndOnes, d$valData[["predTweedie"]])
  })
}

#' AUC of each fold
#'
#' @param rocs Output of `rocPerFold()`.
#' @return Numeric vector; `NA` where the fold has no ROC curve.
aucPerFold <- function(rocs) {
  vapply(rocs, function(r) if (is.null(r)) NA_real_ else as.numeric(r$auc), numeric(1))
}

#' Cross-validation folds for runXGBOOST(), stratified on whether each row is positive.
#'
#' Folds built on the raw response put rare positives in one group, so a single fold could hold every
#' positive and train on zeros only (ELF 3.1.2). Stratified on positive / not, the positives are spread
#' across the folds, and every fold trains with some once there are 2 or more.
#'
#' @param response Numeric response.
#' @param nFolds Number of folds.
#' @return List of `nFolds` lists, each with `keepAll` (all row indices) and `keepEval` (the
#'   fold's held-out row indices).
cvFolds <- function(response, nFolds) {
  folds <- caret::createFolds(factor(response > 0), k = nFolds, list = TRUE, returnTrain = FALSE)
  lapply(folds, function(f) list(keepAll = seq_along(response), keepEval = f))
}

#' Stop when a cross-validation fold would train without any positive observation.
#'
#' xgboost trains each fold's model on `keepAll` without the fold's `eval_set` (`keepEval`). If those
#' rows hold no positive, the model predicts NaN for every row (ELF 3.1.2: one ignition in 2002-2022,
#' which landed in one fold's held-out rows), and pROC::roc() then stops with "No control observation".
#'
#' @param response Numeric response.
#' @param folds Output of `cvFolds()`.
#' @param type `"ignition"` or `"escape"`; used in the error message.
#' @return `NULL`, invisibly. Stops if any fold would train without a positive.
stopIfFoldsLackPositives <- function(response, folds, type) {
  noPositive <- vapply(folds, function(f) !any(response[setdiff(f$keepAll, f$keepEval)] > 0), logical(1))
  if (any(noPositive)) {
    stop("Too few ", type, "s to fit ", length(folds), " cross-validation folds: ",
         sum(response > 0), " positive observation(s), and ", sum(noPositive),
         " fold(s) would train without any.", call. = FALSE)
  }
}

#' Message reporting the mean AUC across folds
#'
#' @param aucs Output of `aucPerFold()`.
#' @return Character string.
meanRocMessage <- function(aucs) {
  if (all(is.na(aucs)))
    return("mean roc:  not computed (no validation fold had both outcomes)")
  paste0("mean roc:  ", format(mean(aucs, na.rm = TRUE), digits = 3),
         if (anyNA(aucs)) paste0(" (over ", sum(!is.na(aucs)), " of ", length(aucs),
                                 " folds)") else "")
}

#' Fit a negative binomial GLM of `ignitions` on all other columns
#'
#' Cannot currently be reached; see `buildModel()`.
#'
#' @param dat data.frame with an `ignitions` column.
#' @return `pROC::roc` object of the in-sample predictions against ignition presence, not the
#'   fitted model.
runGLM.NB <- function(dat) {
  system.time(nb <- glm.nb(ignitions ~., data = dat))
  predNB <- predict(nb, dat, type = "response")
  ignZeroAndOnes <- pmin(1, dat$ign)
  (roc_curveNB <- roc(ignZeroAndOnes, predNB))
}

#' Paste with `_`, for cache function names and filenames
#'
#' @param ... Passed to `paste()`.
#' @param sep Separator.
#' @return Character string.
functionNameHelper <- function(..., sep = "_") {
  paste(..., sep = sep)
}

#' Plot the predicted response against each covariate of one group
#'
#' @param df data.table made in `setupPlots()`, with columns `value`, `predictedProb`, `varFac`.
#' @param fuelOrClimate `"Fuel"` or `"Climate"`: the group of covariates to plot.
#' @param labels Named character vector giving the group of each covariate.
#' @param value Unused; `value` in the plot is the column of `df`.
#' @param jitter Jitter width on x; the height is `jitter / 7000`.
#' @param colors Named vector of colours, one per covariate.
#' @param fuelOrClimateInd Index into `colors` and `labels` of the covariates in this group.
#' @param igOrEsc `"ignition"` or `"escape"`; used in the y-axis label.
#' @return A ggplot.
plotPredictions <- function(df, fuelOrClimate, labels, value, jitter, colors, fuelOrClimateInd, igOrEsc) {

  ggplot(df[varFac %in% names(labels)[labels %in% fuelOrClimate]],
         aes(x = value, y = predictedProb, group = varFac, col = varFac)) +
    geom_point() +
    geom_jitter(width = jitter, height = jitter/7e3) +
    geom_smooth(span = 1) +
    scale_color_manual(aesthetics = "colour", values = colors[fuelOrClimateInd],
                       labels = names(labels)[fuelOrClimateInd],
                       name = "Climate:",
                       guide = guide_legend(reverse = TRUE, title.position = "top", order = 0)) +
    ggplot2::xlab(paste0("Scaled, centred")) +
    ggplot2::ylab(paste0("Predicted probability of ", igOrEsc)) +
    theme_bw()
}

#' Build the response plots for the fitted xgboost models
#'
#' For each fold and covariate, predicts along the covariate's range (limited to -2 to 2, in
#' steps of 0.1, in standardised units) with every other covariate at 0, its mean.
#'
#' @param modelOnly List of the per-fold xgboost models.
#' @param dat data.table of standardised covariates used for the fit.
#' @param igOrEsc `"ignition"` or `"escape"`.
#' @return List of two ggplots, `Climate` and `Fuel`. Covariates whose names match
#'   `CMD|light|positiveCG` are the climate ones.
setupPlots <- function(modelOnly, dat, igOrEsc) {

  cnNoIgnNoEsc <- colnames(dat) |> setdiff(c(fireSenseUtils::ignitionsTxt, fireSenseUtils::escapesTxt, "year", "pixelID"))
  dat <- dat[, ..cnNoIgnNoEsc]

  df <- Map(fold = seq_along(modelOnly), function(fold) {
    df <- Map(nam = cnNoIgnNoEsc, function(nam) {
      rr <- range(dat[, ..nam], na.rm = TRUE) * 10
      rr[1] <- floor(rr[1])
      rr[2] <- ceiling(rr[2])
      df1 <- data.frame(pmax(-20,pmin(20,(rr[1]:rr[2])))/10) |> setNames(nam)
      cnHere <- setdiff(cnNoIgnNoEsc, nam)
      df0 <- lapply(cnHere, function(x) list(0)) |> data.frame() |> setNames(cnHere)
      df <- data.frame(df1, df0)
      # Predict along the regularl sequence of x-axis values
      predictedProb <- predict(modelOnly[[fold]], newdata = df)
      data.table(df1, predictedProb) |>
        melt(measure.vars = nam, stringsAsFactor = FALSE)
    })
    rbindlist(df, use.names = TRUE)
  })
  df <- rbindlist(df, use.names = TRUE)
  set(df, NULL, "variable", as.character(df$variable))
  setorderv(df, "variable")

  set(df, NULL, "varFac", factor(df$variable, levels = cnNoIgnNoEsc)) # keeps colours constant
  set(df, NULL, "varInt", as.integer(df$varFac))

  levels(df$varFac) # gives order that ggplot2 will use
  setorderv(df, "varFac")

  vals <- cnNoIgnNoEsc # stays constant colour regardless of importance
  colors <- NULL
  colOptions <- c("Paired", "Gr")
  colors <- RColorBrewer::brewer.pal(length(vals), "Paired")
  if (length(vals) != length(colors))
    colors <- colorRampPalette(colors)(length(vals))
  names(colors) <- vals
  labels <- rep("Fuel", length(vals))
  names(labels) <- vals
  climateGrep <- "CMD|light|positiveCG"
  climateInd <- grep(climateGrep, vals) # this is how I identify climate vars: not robust!!!!!
  set(df, NULL, "FuelOrClimate", "Fuel")
  set(df, which(df$variable %in% names(labels[climateInd])), "FuelOrClimate", "Climate")

  labels[climateInd] <- "Climate"
  FuelInd <- -climateInd

  jitter <- 0.05

  aa <- Map(foc = c("Climate", "Fuel"), focInd = list(climateInd, FuelInd), function(foc, focInd) {
    plotPredictions(df, fuelOrClimate = foc, labels, value,
                    jitter, colors, fuelOrClimateInd = focInd, igOrEsc)
  })
  aa
}

