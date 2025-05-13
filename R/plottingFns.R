#plotFnLogitIgnition renamed plotFnLogIgnition as logit is incorrect

IgEscapePlots <- function(
    dt = ignitionData$covariates, bestModel = ignitionModel,
    climVar = sim$climateVariablesForFire$ignition,
    fsProcess = "ignition", #or escape
    family = P(sim)$ignitionFamily,
    rescalers = ignitionData$ignitionRescalers, #will be identical to escape
    plotBiomass = P(sim)$plot_fuelBiomassPerPrediction,
    ignitionFitRTM = sim$ignitionFitRTM,
    oPath = outputPath(sim), studyAreaName = P(sim)$.studyAreaName) {

  #general things
  dt <- copy(dt)
  ff <- as.character(bestModel$call$formula)
  y <- bestModel$call$formula[[2L]]
  if (any(c("year", "yr") %in% tolower(names(dt)))) {
    xvar <- intersect(c("year", "yr"), tolower(names(dt)))
  } else {
    xvar <- rows
  }
  formForNull <- as.formula(paste0(ff[[2]], ff[[1]], "1"), env = .GlobalEnv)
  nullModel <- glmmTMB(formForNull, dat = dt, family = eval(family))

  ## https://stackoverflow.com/a/68684973 -- NOT VERY APPROPRIATE FOR RE model
  pseudoR2 <- as.numeric(1 - logLik(bestModel) / logLik(nullModel))

  #plotting
  #TODO: NULL model already produced in buildModel - avoid duplicating it?
  if (terms(nullModel) != terms(bestModel)) {
    #### Plotting glmmTMB ####
    ## build two prediction datasets
    ## both predict across quantiles of climate variable
    ## the first dataset will have class means for cover and forest biomass
    ## the second will have cover values of 1 (i.e. representing complete cover for non-forest)
    ## and biomass values representing the mean for forested pixels (i.e. complete cover of forest)
    ## the 2nd dataset is different because the first one, while more representative of the actual landscape,
    ## implicitly includes non-forest pixels in the calculation of the mean.

    N <- 30
    p <- as.data.table(lapply(dt, function(pp) if (is.numeric(pp)) mean(pp) else pp[1:N]))
    terms <- terms(bestModel)
    termsNonClimate <- setdiff(attr(terms, "term.labels"), climVar)

    ## populate a prediction dataset with quantiles of climate variable and mutually exclusive veg.

    interpolateClimVar <- seq(quantile(dt[[climVar]], 0.1),
                              (quantile(dt[[climVar]], 0.95) * 1.5), length.out = N)
    set(p, NULL, climVar, interpolateClimVar)


    termsNoInteraction <- termsNonClimate[termsNonClimate %in% names(dt)]
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

    termsUsingCover <- as.vector(dt[, lapply(.SD, max), .SDcol = termsNoInteraction])
    termsUsingBiomass <- names(termsUsingCover[termsUsingCover > 1])
    termsUsingCover <- setdiff(names(termsUsingCover), termsUsingBiomass)


    for (val1 in c(termsUsingCover)) {
      set(pAll, which(!pAll$val %in% val1), val1, 0)
    }

    #set minimum biomass as whatever is in data (likely log(100)-1)
    minBiomass <- min(dt[, .SD, .SDcol = termsUsingBiomass])

    for (val1 in c(termsUsingBiomass)) {
      set(pAll, which(!pAll$val %in% c(val1, termsUsingCover)), val1, minBiomass)
    }

    #copy pAll for plot #2 before the data are modified for plot #1
    pAll2 <- copy(pAll)

    #TODO: caching preds is not currently working with reproducible 2.1.2 or 2.1.2.9007 (recursion error)
    system.time({
      preds <- predict(object = bestModel, newdata = pAll, se.fit = TRUE, re.form = NA)
    })
    pAll[, pred := expit(preds$fit)]
    pAll[, val1 := factor(val)]
    pAll[, upper := expit(preds$fit + preds$se.fit)]
    pAll[, lower := expit(preds$fit - preds$se.fit)]

    resInKm2 <- prod(res(ignitionFitRTM)) / 1e6 ## 1e6 m^2 == 1 km^2
    labelToUse <- paste("Ignition rate per", resInKm2, "km^2")
    filenameToUse <- paste0("IgnitionRatePer", resInKm2, "km2_",
                            studyAreaName, "_meanByClass_", climVar)

    titl <- paste0("fireSense_IgnitionFit:", studyAreaName,
                   " (", basename(oPath), ")",
                   " -- Pseudo ")
    titl2 <- paste0(round(pseudoR2, 3))

    if (!is.null(rescalers)) {
      pAll <- rescaleVarsByMagnitude(pAll, 1/rescalers) # invert it
      dt <- rescaleVarsByMagnitude(dt, 1/rescalers) #in case climate is rescaled
    }



    Plots(data = pAll, fn = plotFnLogIgnition, # xColName = colName,
          ggylab = labelToUse,
          subtitle = "using mean cover and biomass per pixel",
          fillTitle = "veg. covariate",
          .plotInitialTime = NULL, # this means "ignore what `.plotInitialTime says; use only .plots`
          # centred = centred,
          climateVar = climVar,
          # origXmax = max(dt[[colName]]), ## if supplied, adds bar to plot
          ggTitle = bquote(.(titl)~R^2 == .(titl2)),
          rawClimate =  dt[[climVar]],
          filename = filenameToUse)

    ## make second prediction using mean forest or alternatively 100% non-forest cover
    if (!is.null(attributes(ignitionFitRTM)$meanForestB) ||
        !is.null(plotBiomass)) {

      Bunit <- ifelse(!is.null(plotBiomass),
                      plotBiomass,
                      log(attributes(ignitionFitRTM)$meanForestB))
      BunitForLabel <- round(exp(Bunit), digits = 0)

      for (val2 in termsUsingBiomass) {
        set(pAll2, which(pAll2$val %in% val2), val2, Bunit)
      }

      for (val2 in termsUsingCover) {
        set(pAll2, which(pAll2$val %in% val2), val2, 1) #set the variable to 1 representing complete cover
      }

      #TODO: caching preds is not currently working with reproducible 2.1.2 or 2.1.2.9007 (recursion error)
      system.time({
        preds <- predict(object= bestModel, newdata = pAll2, se.fit = TRUE, re.form = NA)# |>
        # Cache(omitArgs = "object", .cacheExtra = forms[whBest])
      })

      pAll2[, pred := expit(preds$fit)]
      pAll2[, val1 := factor(val)]
      pAll2[, upper := expit(preds$fit + preds$se.fit)]
      pAll2[, lower := expit(preds$fit - preds$se.fit)]

      if (!is.null(rescalers)) {
        pAll2 <- rescaleVarsByMagnitude(pAll2, 1/rescalers) # invert it
      }

      filenameToUse <- paste0("IgnitionRatePer", resInKm2, "km2_", studyAreaName, "_fullCoverAndBiomass_", climVar)
      Plots(data = pAll2, fn = plotFnLogIgnition,
            ggylab = labelToUse,
            subtitle = paste0("per ", BunitForLabel, " g B/m2 or 100% cover"),
            fillTitle = "veg. covariate",
            .plotInitialTime = NULL, # this means "ignore what `.plotInitialTime says; use only .plots`
            climateVar = climVar,
            rawClimate = dt[[climVar]],
            # origXmax = max(sim$fireSense_ignitionCovariates[[colName]]), ## if supplied, adds bar to plot
            ggTitle = bquote(.(titl)~R^2 == .(titl2)),
            filename = filenameToUse)
    }

    #rescale M once again for this final prediction
    dt <- rescaleVarsByMagnitude(dt, rescalers)

    #TODO: caching preds is not currently working with reproducible 2.1.2 or 2.1.2.9007 (recursion error)
    system.time({
      fittedNoRE <- predict(object = bestModel, newdata = dt, se.fit = FALSE, re.form = NA,
                            type = "response") #|>
      # Cache(.functionName = "predict_forFitted_v_Obs_Ignitions",
      #       omitArgs = "object", .cacheExtra = forms[whBest])
    })

    plotData <- data.table(dt)
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
          ggTitle = paste(studyAreaName, "fireSense_IgnitionFit: obs. vs. fit"),
          ggSubtitle = paste0("Correlation = ", round(correl, 2)),
          filename = paste0(fsProcess, "_NumFiresFitted_", studyAreaName))

  } else {
    message("null ", fsProcess, " model was superior therefore no plotting will occur - please review formula")
  }
}

plotFnLogIgnition <- function(pAll, subtitle = NULL, ggylab,
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

fittedVsObservedPlot <- function(d, ggTitle, ggSubtitle = NULL, ggylab, xColName)  {
  ggplot <- ggplot(data = d, aes_string(x = xColName, y = "value", colour = "variable")) +
    stat_summary(aes(fill = variable), fun.data = mean_ci,
                 geom = "ribbon", alpha = 0.3, show.legend = FALSE) +
    stat_summary(fun = mean, geom = "line", size = 1) +
    scale_color_discrete(labels = c("obsFires" = "observed no. fires",
                                    "predFires" = "fitted no. fires")) +
    theme_bw() +
    theme(legend.position = "bottom") +
    labs(y = ggylab, x = xColName, title = ggTitle,
         subtitle = ggSubtitle, colour = "")
  ggplot
}
