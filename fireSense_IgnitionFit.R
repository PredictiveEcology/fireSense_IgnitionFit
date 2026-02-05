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
  version = list(fireSense_IgnitionFit = "1.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_IgnitionFit.Rmd"),
  loadOrder = list(after = "fireSense_dataPrepFit",
                   before = "fireSense_dataPrepPredict"),
  reqdPkgs = list("data.table", "dplyr", "PredictiveEcology/SpaDES.core@box (>= 2.1.8.9006)",
                  "PredictiveEcology/fireSenseUtils@development (>= 0.0.5.9090)",
                  "glmmTMB", "mirai",
                  "ggplot2", "ggpubr", "MASS", "magrittr",
                  "numDeriv", "parallel", "parallelly",
                  "PredictiveEcology/pemisc@development",
                  "PredictiveEcology/reproducible@AI (>= 2.1.2.9067)",
                  #TODO correct this when reproducible is merged - it is due to cache(predict)
                  "RhpcBLASctl", # "Matrix", # "ModelOriented/EIX",
                  "caret", "pROC",
                  "PredictiveEcology/SHAPforxgboost (>= 0.1.3.9001)", "xgboost (>=3.0.0)", "lightgbm", # install.packages('xgboost', repos = c('https://dmlc.r-universe.dev', 'https://cloud.r-project.org'))
                  # "PredictiveEcology/SpaDES.core@development (>= 2.0.2.9006)",
                  "terra"),
  parameters = bindrows(
    defineParameter("crossValType", "character", c("time-ordered", "crossValidation"), NA, NA,
                    "How the cross validation should happen, time-ordered or regular k-fold crossValidation"),
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
    defineParameter("modelAlgorithm", "character", "xgboost", NA, NA,
                    "Can be `xgboost`, `glmmtmb`, `glm.nb`, `glmmadaptive`, `glm`; only `xgboost` is supported currently"),
    # defineParameter("useMirai", "logical", FALSE, NA, NA,
    #                 "if `TRUE`, then this module will use parallism for fitting and plotting"),
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

frequencyFitRun <- function(sim) {

  #TODO: make cache smart by digest args in advance
  whichProcessesToFit <- P(sim)$whichProcessesToFit
  names(whichProcessesToFit) <- igOrEscNames(whichProcessesToFit, post = "Fitted", case = "sen")

  fireSense_FittedModels <-
    purrr::pmap(.l  = list(igOrEsc = whichProcessesToFit), sim = sim,
                .f = buildModelsFitModels)
  # Map(igOrEsc = whichProcessesToFit, MoreArgs = list(sim = sim),
  #     buildModelsFitModels)

  # put to sim sim$fireSense_IgnitionFitted, sim$fireSense_EscapeFitted
  list2env(fireSense_FittedModels, envir = envir(sim))
  on.exit() # let the mirai keep going, so remove the on.exit close of the mirai object
  return(invisible(sim))

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

prepareCovariates <- function(formula, covariates, rescaleVars, modelAlgorithm) {

  covariates <- copy(setDT(covariates))
  if (any(c("year", "yr") %in% tolower(names(covariates)))) {
    xvar <- intersect(c("year", "yr"), tolower(names(covariates)))
  } else {
    xvar <- rows #TODO what is this?
  }


  if (grepl("xgb", modelAlgorithm) %in% FALSE) {
    formula <- as.formula(formula, env = .GlobalEnv)
    terms <- terms.formula(formula)

    if (attr(terms, "response")) {
      y <- formula[[2L]]
    } else {
      stop("Incomplete formula, the LHS is missing.")
    }
  }

  # if (!is.data.table(covariates))
  #   covariates <- as.data.table(covariates)

  if (rescaleVars) {
    if (grepl("xgb", modelAlgorithm) %in% FALSE) {
      # rescalers <- abs(sapply(covariates[, .SD, .SDcol = toRescale], FUN = max))
      message("Variables outside of [0,10] range will be rescaled to [0,10]")
      toRescale <- setdiff(names(covariates),
                           c("pixelID", ignitionsTxt, escapesTxt, "year", "yearChar"))
      rescalers <- sapply(covariates[, .SD, .SDcol = toRescale], max)
      needRescale <- sapply(rescalers, FUN = function(x) !inRange(x, 0, 10))
      cols <- names(rescalers)[which(needRescale)]
      message("rescaling the following variables: ", paste(cols, collapse = ", "))
      ignitionRescalers <- 10^(floor(log10(abs(rescalers[cols])))) # if range is 0,1, need + 1
      covariates <- rescaleVarsByMagnitude(covariates, ignitionRescalers)
    } else {

      SDcols <- setdiff(colnames(covariates), c("yearChar", ignitionsTxt, "escapes"))
      scaledData <- scale(covariates[, ..SDcols])
      origIgnitions <- covariates[[ignitionsTxt]]
      origEscapes <- covariates[["escapes"]]
      centeringData <- attributes(scaledData)
      covariates <- as.data.table(scaledData)
      set(covariates, NULL, ignitionsTxt, origIgnitions)
      if (!is.null(origEscapes))
        set(covariates, NULL, "escapes", origEscapes)

      # covariates <- covariates[,
      #                          append(
      #                            list(# yearChar = yearChar,
      #                              ignitions = ignitions
      #                              #, nfLCC_100 = nfLCC_100,
      #                              #, nfLCC_50_80 = nfLCC_50_80,
      #                              #, youngAge = youngAge
      #                            ),
      #                            lapply(.SD, scale, scale = TRUE)),
      #                          .SDcols = setdiff(colnames(covariates), c("yearChar", ignitionsTxt))
      #                          #.SDcols = c("CMDsm", "Betu_pap", "Pc_gl.Lr_la", "Pice_mar",
      #                          #"Pn_co.Pn_ba", "Pp_ba.Pp_tr", "lightning")
      # ]
      # cols <- grep("V1", value = TRUE, colnames(covariates))
      # setnames(covariates, old = cols, new = gsub(".V1", "", cols))
      # if (escapesTxt %in% colnames(covariates)) {
      #   set(covariates, NULL, escapesTxt, covariates[[escapesTxt]])
      # }

      setattr(covariates, name = "scaleData", value = centeringData)
      ignitionRescalers <- NULL #so that fire fireSense_IgnitionFit can add it

    }
  } else {
    ignitionRescalers <- NULL #so that fire fireSense_IgnitionFit can add it
  }

  return(list(covariates = covariates,
              formula = formula,
              ignitionRescalers = ignitionRescalers,
              xvar = xvar))
}

buildModel <- function(covariates, formula, type = "ignition",
                       climVar,# = sim$climateVariablesForFire$ignition,
                       family, # useMirai = FALSE,
                       digestOfData,
                       modelAlgorithm, nFolds = 5, crossValType = c("time-ordered", "crossValidation"),
                       formsToRun = c("full", "InterceptOnly", "climateOnly")) {

  # convert to data.table --> easier to work with
  # dat <- as.data.table(covariates)
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
  varToInt <- ignitionsTxt
  if (type == "escape") {
    varToInt <- c(ignitionsTxt, escapesTxt)
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
                              type = type, nFolds = nFolds, crossValType = crossValType[1]) #(nam, form, dat, family, type, climVar)
  } else {
    forms <- forms[formsToRun]

    for (i in c("yearChar"))
      set(dat, NULL, i, factor(dat[[i]]))
    for (ind in seq_along(forms)) {
      nam <- names(forms)[[ind]]
      form <- forms[[ind]]
      if (any(grepl("tmb", modelAlgorithm))) {
        mods[[nam]] <- runGlmmTMB(nam, form, dat, family, type, climVar)
      } else if (any(grepl("adaptive", modelAlgorithm))) {
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



runGlmmTMB <- function(nam, form, dat, family, type, climVar) {

  system.time(po <- glm(ignitions ~ ., data = dat3Forxgboost[-trainIndexK[[1]]], family = "poisson"))
  predGLM <- predict(po, dat3Forxgboost[trainIndexK[[1]], ], type = "response")
  ignZeroAndOnes <- pmin(1, dat3Forxgboost$ign[trainIndexK[[1]]])
  (roc_curveNB <- roc(ignZeroAndOnes, predGLM))

  formSimple <- paste("ignitions ~ ", paste(sep = " + ", "youngAge:CMDsm", "nfLCC_100:CMDsm",
                                            "nfLCC_50_80:CMDsm", "Betu_pap:CMDsm", "Pc_gl.Lr_la:CMDsm",
                                            "Pice_mar:CMDsm", "Pn_co.Pn_ba:CMDsm", "Pp_ba.Pp_tr:CMDsm",
                                            "lightning:CMDsm"))
  formSimple <- as.formula(formSimple, env = .GlobalEnv)
  system.time(po2 <- glm(formSimple, data = dat3Forxgboost[-trainIndexK[[1]]], family = poisson(link = "logit")))
  predGLM2 <- predict(po2, dat3Forxgboost[trainIndexK[[1]], ], type = "response")
  ignZeroAndOnes <- pmin(1, dat3Forxgboost$ign[trainIndexK[[1]]])
  (roc_curveNB <- roc(ignZeroAndOnes, predGLM2))


  predGA <- predict(out[[1]], dat3Forxgboost[-trainIndexK[[1]], ], type = "mean_subject")
  ignZeroAndOnes <- pmin(1, dat3Forxgboost$ign[-trainIndexK[[1]]])
  (roc_curveGA <- roc(ignZeroAndOnes, predGA))

  return(out[[i]])

  print(st2 <- system.time(mu2 <- predict(fit2, newdata = dat, se.fit = FALSE, type_pred = "response", type = "mean_subject")))
  disp2 <- exp(fit2$phis)
  #aa3 <- rnbinom(n = length(mu2), mu = mu2, size = disp2)
  #table(aa3)
  aa2 <- list()
  N <- 100
  for (i in 1:N) aa2[[i]] <- rnbinom(n = length(mu2), mu = mu2, size = log(disp2))
  aa3 <- do.call(cbind, aa2)
  bb3 <- table(aa3)/N
  # cc <- apply(aa1, 2, table)
  sum(as.numeric(names(bb3)) * bb3)



  print(system.time(fit <- glmm.zinb(fixed = ignitions ~  youngAge:CMDsm + nfLCC_100:CMDsm +
                                       nfLCC_50_80:CMDsm + Betu_pap:CMDsm + Pc_gl.Lr_la:CMDsm +
                                       Pice_mar:CMDsm + Pn_co.Pn_ba:CMDsm + Pp_ba.Pp_tr:CMDsm,
                                     random = ~ 1 | yearChar, data = dat, zi_fixed = ~CMDsm, niter  = 100) ))
  #system.time(out <- glmmTMB(form, data = dat,
  #                           ziformula = ziform,
  #                           family = nbinom1(link = "logit")))
  mu <- predict(fit, newdata = dat, level = 0, se.fit = FALSE, type = "response")
  mu <- exp(mu)
  # disp <- sigma(fit)
  disp <- fit$theta
  aa2 <- rnbinom(n = length(mu), mu = mu, size = disp)
  bb2 <- table(aa2)
  sum(as.numeric(names(bb2)) * bb2)
  ccc <- table(dat$ignitions)
  print(sum(as.numeric(names(ccc))*ccc))
  aa <- list()
  N <- 10
  for (i in 1:N) aa[[i]] <- rnbinom(n = length(mu), mu = mu, size = disp)
  aa1 <- do.call(cbind, aa)
  bb <- table(aa1)/N
  cc <- apply(aa1, 2, table)
  sum(as.numeric(names(bb)) * bb)


  print(st4 <- system.time(out <- glmmTMB(form, data = dat,
                                          ziformula = ziform,
                                          family = nbinom1(link = "log"))))
  print(st5 <- system.time(mu1 <- predict(out, newdata = dat, se.fit = FALSE, type = "response")))
  # mu1 <- expit(mu1)
  disp1 <- sigma(out)
  aa3 <- rnbinom(n = length(mu1), mu = mu1, size = 1/disp1)
  aa3t <- table(aa3)
  sum(as.numeric(names(aa3t)) * aa3t)

  # https://win-vector.com/2014/05/30/trimming-the-fat-from-glm-models-in-r/
  out <- local({
    glmmTMB(form, data = dat,
            ziformula = ziform,
            family = eval(family))} |>
      # trimModelObjectForPrediction(origDat = dat)}  |>
      Cache(.functionName = paste0("glmmTMB_for", "_", nam),
            omitArgs = c(formalArgs(glmmTMB), "level", "x", "origDat"),
            # omitArgs = formalArgs(local),
            .cacheExtra = dig),
    # envir = en) # en
    envir = enObjs) # en
  # envir = .GlobalEnv) # en
  ## Use .cacheExtra: there are lots of arguments to glmmTMB that seemed to be "always different"
  #TODO: will nam be an issue if it is identical for escape and ignition models?

  # out2 <- trimModelObjectForPrediction(
  #   out, origDat = dat,
  #   filename = paste0("trimModelObject_", class(out), ".txt")) # |> Cache()

  # a <- identifyEnvs(out, en)
  out
}



buildModelsFitModels <- function(igOrEsc, sim) {
  # function(igOrEsc) {

  covariatesHere <- igOrEscNames(igOrEsc, post = "Covariates") # paste0("fireSense_", igOrEsc, "Covariates")
  objsNeeded <- c(covariatesHere)
  if (grepl("xgb", Par$modelAlgorithm) %in% FALSE) {
    formulaHere <- igOrEscNames(igOrEsc, post = "Formula") # paste0("fireSense_", igOrEsc, "Formula")
    objsNeeded <- c(objsNeeded, formulaHere)
    modelFamily <- igOrEscNames(igOrEsc, pre = "", post = "Family") # paste0(igOrEsc, "Family")
    family <- P(sim)[[modelFamily]]
    formulaHere <- sim[[formulaHere]] # pull from simList
  } else {
    formulaHere <- NULL
    family <- NULL
  }

  # This takes time for large datasets
  digestOfData <- .robustDigest(mget(objsNeeded, envir = envir(sim)))

  data <- prepareCovariates(formula = formulaHere,
                            covariates = sim[[covariatesHere]],
                            rescaleVars = P(sim)$rescaleVars,
                            modelAlgorithm = Par$modelAlgorithm) |>
    Cache(omitArgs = c("covariates", "formula"), .cacheExtra = digestOfData)


  if (identical(igOrEsc, "escape")) {
    # Escape should not have lightning
    # if don't explicitly copy, then Cache above returns the "lightning"-removed data.table
    data$covariates <- data.table::copy(data$covariates)
    set(data$covariates, NULL, "lightning", NULL)
  }
  nFolds <- 5

  crossValType <- Par$crossValType

  # opt <- options(reproducible.cacheSaveFormat = "rds")
  # on.exit(options(opt)) # redundant; but necessary if it fails during fit
  modelHere <- buildModel(covariates = data$covariates,
                          climVar = sim$climateVariablesForFire$ignition,
                          formula= data$formula,
                          type = igOrEsc, nFolds = nFolds,
                          # type = "ignition",
                          crossValType = crossValType,
                          family = family, # useMirai = P(sim)$useMirai,
                          digestOfData = digestOfData,
                          modelAlgorithm = Par$modelAlgorithm,
                          formsToRun = c("full", "InterceptOnly", "climateOnly")[1]) |>
    Cache(omitArgs = c("covariates", "formula"),
          .functionName = paste0("buildModel for ", igOrEsc),
          .cacheExtra = c(digestOfData, 1), cacheSaveFormat = "rds") # doesn't work with qs
  # options(opt) # redundant; but necessary so stuff below has qs (or original csf)
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
    # library(mirai)
    # This says, "create 2 workers", no more. So, no matter how many mirai are started, they
    #   just stay at 2 cores max.
    # on.exit(daemons(0), add = TRUE)
    # try(daemons(2, dispatcher = FALSE), silent = TRUE) # this is ignored the 2nd time
    message("Plotting ", igOrEsc, "...")
    # argsForDigest <- list(climVar = sim$climateVariablesForFire[[igOrEsc]],
    #                       # rescalers = data$ignitionRescalers,
    #                       fsProcess = igOrEsc, family = P(sim)[[modelFamily]],
    #                       plotBiomass = P(sim)$plot_fuelBiomassPerPrediction,
    #                       ignitionFitRTM = sim$ignitionFitRTM,
    #                       studyAreaName = P(sim)$.studyAreaName,
    #                       oPath = outputPath(sim), IgEscapePlots = IgEscapePlots,
    #                       digestOfData = digestOfData)
    # digestForPlot <- .robustDigest(argsForDigest)
    # argsAll <- append(list(data = data, bestModel = modelHere, digestForPlot = digestForPlot),
    #                   argsForDigest)

    # mir[[igOrEsc]] <- mirai(.expr = {
    #   library(glmmTMB)
    #   library(reproducible)
    #   IgEscapePlots(dt = data$covariates, bestModel = bestModel,
    #                 climVar = climVar,
    #                 rescalers = data$ignitionRescalers,
    #                 fsProcess = fsProcess, family = family,
    #                 plotBiomass = plotBiomass,
    #                 ignitionFitRTM = ignitionFitRTM,
    #                 studyAreaName = studyAreaName,
    #                 oPath = oPath) |>
    #     Cache(omitArgs = formalArgs(IgEscapePlots), .cacheExtra = append(digestForPlot, digestOfData))
    # },
    # .args = argsAll)  # this tells mirai which objects are needed

    # cnNoIgnNoEsc <- colnames(data$covariates) |> setdiff(c(ignitionsTxt, escapesTxt, "year", "pixelID"))
    # dat <- data$covariates[, ..cnNoIgnNoEsc]

    figPath <- figurePath(sim)
    digModels <- attr(modelHere, "tags")# .robustDigest(modelHere)
    modelOnly <- modelHere[grep("Fold", names(modelHere))]
    aa <- setupPlots(modelOnly, dat = data$covariates, igOrEsc) |>
      Cache(omitArgs = c("dat", "modelOnly"),
            .cacheExtra = list(digestOfData = digestOfData, digModels = digModels, plotPredictions = plotPredictions),
            .functionName = paste0(".functionName_", igOrEsc))
    fn1 <- functionNameHelper("FuelClimate", ifelse(igOrEsc == "escape", "", "Lightning"), "predicted", igOrEsc, Par$crossValType[1])
    fn <- functionNameHelper(fn1, format(Sys.time()))
    bb <- quote(ggarrange(plotlist = aa))

    a <- Plots(bb, path = figPath,
               filename = fn,
               type = "png", ggsaveArgs = list(width = 8, height = 5, scale = 1.5)) |>
      Cache(omitArgs = c("data", "filename"),
            .cacheExtra = list(digestOfData = digestOfData, digModels = digModels, filename = fn1),
            .functionName = paste0(".functionName_", igOrEsc)) |> reproducible:::suppressWarningsSpecific("appears to have a much larger size on disk than in memory")

    if (FALSE)
      IgEscapePlots(dt = data$covariates, bestModel = modelOnly,
                    modelAlgorithm = Par$modelAlgorithm,
                    climVar = sim$climateVariablesForFire[[igOrEsc]],
                    rescalers = data$ignitionRescalers,
                    fsProcess = "ignition", family = P(sim)$ignitionFamily,
                    plotBiomass = P(sim)$plot_fuelBiomassPerPrediction,
                    ignitionFitRTM = sim$ignitionFitRTM,
                    studyAreaName = P(sim)$.studyAreaName,
                    oPath = outputPath(sim)) |>
      Cache(omitArgs = c("dt", "bestModel"), .cacheExtra = digestOfData)
  }
  list(modelList = modelList, scaleData = attr(data$covariates, "scaleData"))
  #  }
}

trimModelObjectForPrediction <- function(x, origDat, filename) {
  opts <- options(reproducible.showSimilar = FALSE)
  on.exit(options(opts))
  headDat <- head(origDat)
  dig <- .robustDigest(headDat)

  rr <- try(ranef(x))
  frame <- x$frame
  if (is(rr, "try-error")) {
    keep <- sample(NROW(x$frame), size = 20)
  } else {
    # keep 1 row with each ranef
    keep <- as.data.table(frame)[, .I[1], by = names(rr$cond)]$V1
  }
  code <- character()
  if (!file.exists(filename)) {

    if ( (is.environment(x) || is.list(x) ) ) {
      for (y in setdiff(names(x), "")) {
        if ((is(x[[y]], "list") || is.call(x[[y]]) || is.environment(x[[y]])) ) {
          for (z in setdiff(names(x[[y]]), "")) {
            if ((is(x[[y]][[z]], "list") || is.call(x[[y]][[z]]) || is.environment(x[[y]][[z]])) ) {
              for (w in setdiff(names(x[[y]][[z]]), "")) {
                orig <- x[[y]][[z]][[w]]
                x[[y]][[z]][[w]] <- c()
                pre <- try(predict(x, newdata = headDat) |>
                             Cache(omitArgs = c("newdata"), .cacheExtra = dig))
                theSymb <- paste0(y, "$", z, "$", w)
                if (is(pre, "try-error")) {
                  message("Failed: --------------> ", theSymb)
                  x[[y]][[z]][[w]] <- orig
                } else {

                  code <- c(code, paste0("out$", theSymb, " <- list()"))
                  message("Replaced: ", theSymb)
                }
              }
            } else if (is(x[[y]][[z]], "data.frame")) {
              x[[y]][[z]] <- head(x[[y]][[z]])
              pre <- try(predict(x, newdata = headDat) |>
                           Cache(omitArgs = c("newdata"), .cacheExtra = dig))
              theSymb <- paste0(y, "$", z)
              if (is(pre, "try-error")) {
                message("Failed: --------------> ", theSymb)
                x[[y]][[z]] <- orig
              } else {
                code <- c(code, paste0("out$", theSymb, " <- list()"))
                message("Replaced: ", )
              }
            } else {
              orig <- x[[y]][[z]]
              x[[y]][[z]] <- c()
              pre <- try(predict(x, newdata = headDat) |>
                           Cache(omitArgs = c("newdata"), .cacheExtra = dig))
              theSymb <- paste0(y, "$", z)
              if (is(pre, "try-error")) {
                message("Failed: --------------> ", theSymb)
                x[[y]][[z]] <- orig
              } else {
                code <- c(code, paste0("out$", theSymb, " <- list()"))
                message("Replaced: ", theSymb)
              }
            }
          }
        } else if (is(x[[y]], "data.frame")) {
          x[[y]] <- head(x[[y]])
          pre <- try(predict(x, newdata = headDat) |>
                       Cache(omitArgs = c("newdata"), .cacheExtra = dig))
        } else {
          orig <- x[[y]]
          x[[y]] <- c()
          pre <- try(predict(x, newdata = head(dat)))
          theSymb <- paste0(y)
          if (is(pre, "try-error")) {
            message("Failed: --------------> ", theSymb)
            x[[y]] <- orig
          } else {
            code <- c(code, paste0("out$", theSymb, " <- list()"))
            message("Replaced: ", theSymb)
          }
        }
      }
    }
    code <- c("a <- function(out) {", code, "out", "}")
    cat(code, sep = "\n", file = filename)
  } else {
    fn <- eval(parse(filename), envir = .GlobalEnv)
    x2 <- fn(x)
  }

  x$frame <- frame[keep, ]
  x
}

# MapRunGlmmTMB <- function(ind, forms, dat, family, type, climVar) {
#   nam <- names(forms)[[ind]]
#   form <- forms[[ind]]
#   mod <- runGlmmTMB(nam, form, dat, family, type, climVar)
# }


## Below:  These are likely not needed
os <- function(x, level = 2) {
  if (is.list(x) && level > 0) {
    Map(x = x, function(x) os(x, level = level - 1))
  } else if (is.environment(x))  {
    os(as.list(x, all.names = TRUE), level = level - 1)
    # Map(x = x, function(x) os(as.list(x)[-1]))
  } else {
    obj_size(x)
  }
}
os2 <- function(x, level = 2) {
  os(x) |> unlist() |> sort()
}
notFns <- function(x, level = 2) {
  if (is.list(x) && level > 0) {
    out <- Map(y = as.list(x), function(y) if(is.function(y)) NULL else notFns(y, level = level - 1))
  } else if (is.environment(x) && level > 0) {
    out <- notFns(as.list(x), level = level - 1)
  } else {
    out <- os(x, level = 0)
  }
  out <- out[!sapply(out, is.null)]
  out |> unlist() |> sort()
}



runGLMMAdaptiveWithSimplifications <- function(nam, form, dat, family, type, climVar) {
  en <- new.env(parent = asNamespace("fireSenseUtils"))
  enObjs <- new.env(parent = en)
  if (type == "ignition") {
    ziform <- paste0("~", paste0(climVar, collapse = "+"))
  } else {
    ziform <- "~0"
  }
  # objNamesOutside <- c("dat", "family")
  # objNamesInside <- c("form", "ziform", "type", "nam", "dat", "family", "climVar")
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

    # set.seed(123)
    # dat <- dat[sample(NROW(dat), 100000),] # ELIOT TO REMOVE

    # Put the objects and functions in the other environnemt
    list2env(objs, envir = enObjs)

    # message("Running GLMMadaptive::mixed_model with Zero-Inflated, Mixed effect, Negative Binomial, using:\n",
    #         messageFormulaFn(form))
    oo <- capture.output(objs$model)
    message("Running ", oo, "\n", "where 'form' is: ", messageFormulaFn(form))

    if (FALSE) {

      Require::Require("SimonDedman/gbm.auto")
      dat <- readRDS("~/dat2.rds")
      #sam <- sample(NROW(dat), size = 1e5)
      # Require::Require("gbm.auto")
      # dat <- readRDS("c:/Users/emcintir/dat.rds")
      # dat <- datAll
      sam <- seq(NROW(dat))
      #sample(NROW(dat), size = 1e4)
      dir.create("outputs", showWarnings = FALSE)
      print(st <- system.time(aa <- gbm.auto(samples =as.data.frame(dat[sam]),
                                             resvar = 4,
                                             expvar = setdiff(seq(ncol(dat)), c(1,4)),
                                             ZI = TRUE,
                                             fam2 = "poisson",
                                             # multiplot = FALSE, MLEvaluate = FALSE,
                                             gaus = FALSE,
                                             savedir = "outputs")))

    }

    library(GLMMadaptive)
    print(system.time(
      out[[i]] <- try(local(
        eval(model) |>
          Cache(.functionName = paste0("mixed_model_for", "_", nam),
                omitArgs = c(formalArgs(GLMMadaptive::mixed_model), "level", "x", "origDat"),
                # omitArgs = formalArgs(local),
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
      # objs$dig$thisFamily <- .robustDigest(objs$thisFamily)
      objs$dig$model <- .robustDigest(objs$model)
    }

    objs$dig$thisFamily <- .robustDigest(objs$thisFamily)

  }
  out
}




runXGBOOST <- function(dat, dig, type = "ignition", nFolds = 5,
                       crossValType = c("time-ordered", "crossValidation")) {
  dat3Forxgboost <- dat# [sample(NROW(dat), size = 1e6)]
  # tt <- table(datForXGBoost$ignitions); sum(as.numeric(names(tt)) * tt)

  # Add dummy variables for factor columns -- i.e., the random effects
  if (all(sapply(dat3Forxgboost, is.numeric)) %in% FALSE)
    dat3Forxgboost <- #{
      model.matrix(~ . + 0, data = dat3Forxgboost) |>
      Cache(omitArgs = c("object", "data", "x"), .cacheExtra = dig)# Creates dummy variables


  savedSeed <- .Random.seed
  on.exit(assign(".Random.seed", savedSeed, envir = .GlobalEnv), add = TRUE)
  set.seed(12345) # so kfolds are same, so Caching works correctly below; if dat3 changes number of rows,
  # it will be a totally different sequence; but it will be the same sequence
  # if number of rows doesn't change

  # Setup k-folds
  # nFolds <- 5
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
    trainIndexK <- caret::createFolds(dat3Forxgboost[, ignitions], k = nFolds, list = TRUE, returnTrain = FALSE)
    trainIndexK <- Map(tr = trainIndexK, function(tr) {
      list(seq(NROW(dat3Forxgboost)), tr) |> setNames(indexNames)
    })
  }

  colOrder <- setdiff(colnames(dat3Forxgboost), c("pixelID", "year"))
  colOrder <- colOrder[colOrder %in% grep("yearChar", colnames(dat3Forxgboost), invert = TRUE, value = TRUE)]
  # colOrder <- setdiff(colOrder, "CMDsp") # remove CMDsp
  colOrder <- sample(colOrder)
  dat3Forxgboost <- dat3Forxgboost[, ..colOrder]
  dig <- .robustDigest(dat3Forxgboost)
  # sim <- get("sim", whereInStack("sim")) # work around--> correct way is to pass figPath
  # figPath <- figurePath(sim)
  # rm(sim)
  colnamesNoIgn <- grep(paste0(ignitionsTxt,"|",escapesTxt), colnames(dat3Forxgboost), value = TRUE, invert = TRUE) |>
    sort() # make alphabetical
  dat3ForxgboostNoIgn <- dat3Forxgboost[, ..colnamesNoIgn]

  ignOrEscapeColName <- grep(value = TRUE,type, colnames(dat3Forxgboost))
  # if (packageVersion("xgboost") <= "3.0.2.1") {
  #   install.packages('xgboost', repos = c('https://dmlc.r-universe.dev', 'https://cloud.r-project.org'))
  #   stop("Please restart R")
  # }
  st <- system.time(
    mm <- purrr::pmap(
      list(valInd = trainIndexK, kFold = seq(nFolds)),
      function(valInd, kFold)#, indexName = indexNames, digInner = dig, dat3ForxgboostInner = dat3Forxgboost,
        #         dat3ForxgboostNoIgnInner = dat3ForxgboostNoIgn)
      {

        wholeDataset <- valInd[[indexNames[[1]]]]
        valInd <- valInd[[indexNames[[2]]]]
        digValInd <- .robustDigest(valInd) # should be eval set

        # xgboost objects do not save with `qs` ... must be `rds`
        # opt <- options(reproducible.cacheSaveFormat = "rds")
        # on.exit(options(opt)) # redundant; but necessary if it fails during fit
        mTweedie <- xgboost(x = dat3ForxgboostNoIgn[wholeDataset], # should be whole set
                            # dat3Forxgboost[, get(ignOrEscapeColName)][wholeDataset],
                            y = dat3Forxgboost[, get(ignOrEscapeColName)][wholeDataset],
                            # objective = "count:poisson", #
                            objective = "reg:tweedie",
                            nthread = 10,
                            eval_set = valInd, # 0.2, # should be keepEval
                            monitor_training = TRUE,
                            verbosity = 0,
                            # eval_metric = c("auc", "logloss"),
                            nrounds = 100,
                            max_depth = 2,
                            reg_lambda = 0.5,
                            learning_rate = 0.15
        ) |> Cache(omitArgs = c("x", "y", "eval_set"),
                   # ) |> list()} |> Cache(omitArgs = c("x", "y", "eval_set", "..."),
                   .functionName = functionNameHelper("xgboost", type, kFold),
                   .cacheExtra = c(dig, digValInd, type),
                   cacheSaveFormat = "rds")
        # options(opt)
        # mTweedie <- mTweedie[[1]] # remove the list
        # Predict probabilities
        valData <- dat3Forxgboost[valInd, ]
        pred2 <- predict(mTweedie, valData)
        valData <- cbind(valData, predTweedie = pred2)
        if (FALSE) {
          ignZeroAndOnes <- pmin(1, valData[, ignitions])
          # (roc_curvePoisson <- roc(ignZeroAndOnes, d$valData$predPoisson))
          (roc_curveTweedie <- roc(ignZeroAndOnes, valData[["predTweedie"]]))
        }

        # class(mTweedie) <- c("lgb.Booster", class(mTweedie)) # bug in next function; needs to be this class
        shap_values <- shap.values(mTweedie, dat3ForxgboostNoIgn) |>
          Cache(omitArgs = formalArgs(shap.values),
                .functionName = functionNameHelper("shap.values", type, kFold),
                .cacheExtra = c(dig, digValInd, type))
        shapContrib <- shap_values$shap_score
        shapContrib <- shapContrib[, -"(Intercept)"]
        shap_long <- shap.prep(# xgb_model = mTweedie,
          shap_contrib = shapContrib, X_train = dat3ForxgboostNoIgn)  |>
          Cache(omitArgs = formalArgs(shap.prep),
                .functionName = functionNameHelper("shap.prep", type, kFold),
                .cacheExtra = c(dig, digValInd, type))
        list(valData = valData, mod = mTweedie, shap_long = shap_long)
      })
  )
  # stStart <- Sys.time()
  # m <- mirai::mirai(
  #   .args = list(dig = dig, shap_long = shap_long, basefilename = figPath, kFold = kFold,
  #                digValInd = digValInd),
  #   {
  #     library(reproducible)
  #     library(SpaDES.core)
  #     library(SHAPforxgboost)
  #     library(ggplot2)
  #
  #     # These plots take 45 minutes for dataset with 4.5M rows
  #     # plt <- Plots(shap.plot.summary(shap_long),
  #     #              path = basefilename, filename = paste0("Ignitions_SHAP model_", kFold),
  #     #              type = "png") |>
  #     #   Cache(omitArgs = "data", .cacheExtra = c(dig, digValInd),
  #     #         ,
  #     #         .functionName = "shap.plot.summary")
  #
  #     plt2 <- Plots(
  #       shap.plot.dependence(data_long = shap_long,
  #                            x = 'CMDsm',
  #                            y = 'Pice_mar',
  #                            color_feature = 'Column_WV') +
  #         ggtitle("SHAP values of Pice_var vs. CMDsm"),
  #
  #       path = basefilename, filename = paste0("Ignitions_SHAP model_dependencies", kFold),
  #       types = "png") |>
  #       Cache(omitArgs = "data", .cacheExtra = c(dig, digValInd, "plotDep"),
  #             .functionName = "shap.plot.dependence")
  #
  #     # g2 <- shap.plot.dependence(data_long = shap_long, x = 'dayint', y = 'Column_WV', color_feature = 'Column_WV') +  ggtitle("(B) SHAP values of CWV vs. Time trend")
  #
  #   })
  # stEnd <- Sys.time()

  #       # set(valData, NULL, c("predTweedie"), list(pred2))
  #
  #       list(valData = valData, mod = mTweedie)
  #     })
  # )

  # predGlm <- predict(mGlm, valData)
  # cn <- colnames(dat3Forxgboost)
  # # Map(nam = cn, function(nam) data.frame(0) |> setNames(cn))
  #
  # df <- Map(fold = seq(nFolds), function(fold) {
  #   df <- Map(nam = setdiff(cn, c(ignitionsTxt, escapesTxt)), function(nam) {
  #     rr <- range(dat3Forxgboost[, ..nam], na.rm = TRUE) * 10
  #     rr[1] <- floor(rr[1])
  #     rr[2] <- ceiling(rr[2])
  #     df1 <- data.frame(pmax(-20,pmin(20,(rr[1]:rr[2])))/10) |> setNames(nam)
  #     cnHere <- setdiff(cn, nam)
  #     df0 <- lapply(cnHere, function(x) list(0)) |> data.frame() |> setNames(cnHere)
  #     df <- data.frame(df1, df0)
  #     # Predict along the regularl sequence of x-axis values
  #     predictedProb <- predict(mm[[fold]]$mod, newdata = df)
  #     data.table(df1, predictedProb) |>
  #       melt(measure.vars = nam, stringsAsFactor = FALSE)
  #   })
  #   rbindlist(df, use.names = TRUE)
  # })
  # df <- rbindlist(df, use.names = TRUE)
  # set(df, NULL, "variable", as.character(df$variable))
  # setorderv(df, "variable")
  #
  # df2 <- Map(m = mm, function(m) xgb.importance(m$mod)) |> rbindlist()
  # importance <- df2[, lapply(.SD, mean), by = "Feature"] |> setorderv("Gain", order = -1L)
  #
  # set(df, NULL, "ord", match(df$variable, importance$Feature))
  # setorderv(df, "ord", order = -1L)
  if (FALSE)  {
    # ggIgnit <-
    ggplot(df) +
      aes(x = value, y = predictedProb, group = variable, col = variable) +
      # ggIgnit <- ggplot(df) + aes(x = value, y = predictedProb, group = variable, col = variable) +
      geom_smooth(span = 1) +
      # ggplot2::geom_jitter(height = 0.00002, size = 0.5) +
      geom_point(size = 0.3) +
      # geom_line() +
      theme_bw() +
      scale_color_brewer(palette="Paired") +
      ggplot2::xlab(paste0("Scaled, centred")) +
      ggplot2::ylab(paste0("Predicted probability of ", type)) +
      guides(col = guide_legend(reverse = FALSE, col = factor(df$predictedProb)))
  }

  #
  #   set(df, NULL, "varFac", factor(df$variable, levels = colnamesNoIgn)) # keeps colours constant
  #   set(df, NULL, "varInt", as.integer(df$varFac))
  #
  #   levels(df$varFac) # gives order that ggplot2 will use
  #   setorderv(df, "varFac")
  #
  #   vals <- unique(importance$Feature)
  #   vals <- levels(df$varFac) # MUST USE THIS FOR GGPLOT2 TO GET CORRECT LABELS
  #   # vals <- colnamesNoIgn # stays constant colour regardless of importance
  #   colors <- RColorBrewer::brewer.pal(length(vals), "Paired")
  #   names(colors) <- vals
  #   labels <- rep("Fuel", length(vals))
  #   names(labels) <- vals
  #   climateGrep <- "CMD|light"
  #   climateInd <- grep(climateGrep, vals) # this is how I identify climate vars: not robust!!!!!
  #   set(df, NULL, "FuelOrClimate", "Fuel")
  #   set(df, which(df$variable %in% names(labels[climateInd])), "FuelOrClimate", "Climate")
  #
  #   labels[climateInd] <- "Climate"
  #   FuelInd <- -climateInd

  # a <- ggplot(df, aes(x = value, y = predictedProb, group = varFac, col = varFac)) +
  #   geom_point(# data = df[varFac %in% names(labels)[labels %in% "Climate"]],
  #     # aes(x = value, y = predictedProb, group = varFac, col = varFac)
  #   ) +
  #   geom_jitter(width = jitter, height = jitter/7e3) +
  #   geom_smooth(# data = df[varFac %in% names(labels)[labels %in% "Climate"]],
  #     # aes(x = value, y = predictedProb, group = varFac, col = varFac),
  #     span = 1) +
  #   # scale_color_manual(aesthetics = "colour", values = colors[climateInd],
  #   #                    labels = names(labels)[climateInd],
  #   #                    # breaks = names(colors)[climateInd],
  #   #                    name = "Climate:",
  #   #                    guide = guide_legend(reverse = TRUE, title.position = "top", order = 0)) +
  #   ggplot2::xlab(paste0("Scaled, centred")) +
  #   ggplot2::ylab(paste0("Predicted probability of ", type)) +
  #   theme_bw() +
  #   scale_color_brewer(palette="Paired") +
  #   facet_grid(cols = vars(FuelOrClimate))
  #
  # ggplot() +
  #   geom_point(data = df[varFac %in% names(labels)[labels %in% "Climate"]],
  #              aes(x = value, y = predictedProb, group = varFac, col = varFac)) +
  #   geom_smooth(data = df[varFac %in% names(labels)[labels %in% "Climate"]],
  #               aes(x = value, y = predictedProb, group = varFac, col = varFac), span = 1) +
  #   scale_color_manual(aesthetics = "colour", values = colors[climateInd],
  #                      labels = names(labels)[climateInd],
  #                      # breaks = names(colors)[climateInd],
  #                      name = "Climate:",
  #                      guide = guide_legend(reverse = TRUE, title.position = "left", order = 0)) +
  #   new_scale_colour() +
  #   geom_point(data = df[varFac %in% names(labels)[labels %in% "Fuel"]], aes(x = value, y = predictedProb, group = varFac, col = varFac)) +
  #   geom_smooth(data = df[varFac %in% names(labels)[labels %in% "Fuel"]], aes(x = value, y = predictedProb, group = varFac, col = varFac), span = 1) +
  #   scale_color_manual(aesthetics = "colour",
  #                      values = colors[-climateInd],
  #                      labels = names(labels)[-climateInd],
  #                      name = "Fuels:",
  #                      guide = guide_legend(reverse = TRUE, title.position = "left", order = 1)) +
  #   theme_bw()


  # jitter <- 0.05
  # a <- ggplot(df[varFac %in% names(labels)[labels %in% "Climate"]],
  #             aes(x = value, y = predictedProb, group = varFac, col = varFac)) +
  #   geom_point(# data = df[varFac %in% names(labels)[labels %in% "Climate"]],
  #              # aes(x = value, y = predictedProb, group = varFac, col = varFac)
  #              ) +
  #   geom_jitter(width = jitter, height = jitter/7e3) +
  #   geom_smooth(# data = df[varFac %in% names(labels)[labels %in% "Climate"]],
  #               # aes(x = value, y = predictedProb, group = varFac, col = varFac),
  #               span = 1) +
  #   scale_color_manual(aesthetics = "colour", values = colors[climateInd],
  #                      labels = names(labels)[climateInd],
  #                      # breaks = names(colors)[climateInd],
  #                      name = "Climate:",
  #                      guide = guide_legend(reverse = TRUE, title.position = "top", order = 0)) +
  #   ggplot2::xlab(paste0("Scaled, centred")) +
  #   ggplot2::ylab(paste0("Predicted probability of ", type)) +
  #   theme_bw()
  #

  rocs <- lapply(mm, function(d) {
    ignZeroAndOnes <- pmin(1L, d$valData[[ignOrEscapeColName]])
    # (roc_curvePoisson <- roc(ignZeroAndOnes, d$valData$predPoisson))
    (roc_curveTweedie <- pROC::roc(ignZeroAndOnes, d$valData[["predTweedie"]]))
    roc_curveTweedie
  })

  tweedie <- mapply(r = rocs, function(r) {
    as.numeric(r$auc)
  })
  # poiss <- mapply(r = rocs, function(r) {
  #   as.numeric(r$roc_curvePoisson$auc)
  # })
  print(paste("mean roc: ", format(mean(tweedie), digits = 3)))
  # mean(poiss)
  mm2 <- Map(m = mm, function(m) m$mod)
  mm2 <- append(mm2, list(rocs = rocs))

  return(mm2)
}


runGLM.NB <- function(dat) {
  system.time(nb <- glm.nb(ignitions ~., data = dat))
  predNB <- predict(nb, dat, type = "response")
  ignZeroAndOnes <- pmin(1, dat$ign)
  (roc_curveNB <- roc(ignZeroAndOnes, predNB))
}

ignitionsTxt <- "ignitions"
escapesTxt <- "escapes"

functionNameHelper <- function(..., sep = "_") {
  paste(..., sep = sep)
  #   paste(fnName, type, kFold, sep = sep)
}


xgboost2 <- function(...) {
  xgboost(...) |> list()
  # dat3ForxgboostNoIgn[wholeDataset], # should be whole set
  # # dat3Forxgboost[, get(ignOrEscapeColName)][wholeDataset],
  # dat3Forxgboost[, get(ignOrEscapeColName)][wholeDataset],
  # # objective = "count:poisson", #
  # objective = "reg:tweedie",
  # nthread = 10,
  # eval_set = valInd, # 0.2, # should be keepEval
  # monitor_training = TRUE,
  # verbosity = 0,
  # # eval_metric = c("auc", "logloss"),
  # nrounds = 100,
  # max_depth = 2,
  # reg_lambda = 0.5,
  # learning_rate = 0.15
}



plotPredictions <- function(df, fuelOrClimate, labels, value, jitter, colors, fuelOrClimateInd, igOrEsc) {

  ggplot(df[varFac %in% names(labels)[labels %in% fuelOrClimate]],
         aes(x = value, y = predictedProb, group = varFac, col = varFac)) +
    geom_point(# data = df[varFac %in% names(labels)[labels %in% fuelOrClimate]],
      # aes(x = value, y = predictedProb, group = varFac, col = varFac)
    ) +
    geom_jitter(width = jitter, height = jitter/7e3) +
    geom_smooth(# data = df[varFac %in% names(labels)[labels %in% fuelOrClimate]],
      # aes(x = value, y = predictedProb, group = varFac, col = varFac),
      span = 1) +
    scale_color_manual(aesthetics = "colour", values = colors[fuelOrClimateInd],
                       labels = names(labels)[fuelOrClimateInd],
                       # breaks = names(colors)[fuelOrClimateInd],
                       name = "Climate:",
                       guide = guide_legend(reverse = TRUE, title.position = "top", order = 0)) +
    ggplot2::xlab(paste0("Scaled, centred")) +
    ggplot2::ylab(paste0("Predicted probability of ", igOrEsc)) +
    theme_bw()
}



setupPlots <- function(modelOnly, dat, igOrEsc) {

  cnNoIgnNoEsc <- colnames(dat) |> setdiff(c(ignitionsTxt, escapesTxt, "year", "pixelID"))
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

  # vals <- unique(importance$Feature)
  # vals <- levels(df$varFac) # MUST USE THIS FOR GGPLOT2 TO GET CORRECT LABELS
  vals <- cnNoIgnNoEsc # stays constant colour regardless of importance
  colors <- NULL
  colOptions <- c("Paired", "Gr")
  # if(length(vals) != length(colors)) {
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

