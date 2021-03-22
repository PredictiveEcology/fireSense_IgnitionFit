defineModule(sim, list(
  name = "fireSense_IgnitionFit",
  description = paste("Fit statistical models that can be used to parameterize (calibrate)",
                      "the fire ignition component of landscape fire models (e.g. fireSense)."),
  keywords = c("fire frequency", "optimization", "additive property", "poisson",
               "negative binomial", "fireSense"),
  authors = c(
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre")),
    person("Ian", "Eddy", email = "ian.eddy@canada.ca", role = c("aut")),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(),
  version = list(SpaDES.core = "0.1.0", fireSense_IgnitionFit = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_IgnitionFit.Rmd"),
  reqdPkgs = list("DEoptim", "dplyr", "ggplot2", "MASS", "magrittr", "numDeriv", "parallel", "pemisc",
                  "PredictiveEcology/fireSenseUtils@development (>=0.0.4.9048)",
                  "PredictiveEcology/SpaDES.core@development (>=1.0.6.9019)"), # need Plots stuff
  parameters = bindrows(
    defineParameter("autoRefit", c("logical", "character"), default = TRUE, min = NA, max = NA,
                    desc = paste("If the objective function results in a singularity or non-convergence with full ",
                                 "model, should the module cull all effects that are causing the problem and ",
                                 "retry?")),
    defineParameter("cores", "integer", default = 1,
                    desc = paste("non-negative integer. Defines the number of logical cores",
                                 "to be used for parallel computation.",
                                 "The default value is 1, which disables parallel computing.")),
    defineParameter("data", "character", default = "fireSense_ignitionCovariates",
                    desc = paste("a character vector indicating the names of objects in the",
                                 "`simList` environment in which to look for variables present",
                                 "in the model formula. `data` objects should be data.frames.")),
    defineParameter("family", "function, character", default = quote(poisson(link = "identity")),
                    desc = paste("a family function (must be wrapped with `quote()`) or a",
                                 "character string naming a family function.",
                                 "For additional details see `?family`.")),
    defineParameter("fireSense_ignitionFormula", "character", default = NA,
                    desc = paste("formula - as a character - describing the model to be fitted.",
                                 "Piece-wised terms can be specifed using `pw(variableName, knotName)`.")),
    defineParameter("iterDEoptim", "integer", default = 500,
                    desc = "maximum number of iterations allowed (DEoptim optimizer)."),
    defineParameter("iterNlminb", "integer", default = 500,
                    desc = paste("if start is not supplied, iterNlminb defines the number of trials,",
                                 "or searches, to be performed by the nlminb optimizer in order to",
                                 "find the best solution.")),
    defineParameter("lb", "list", default = NULL,
                    desc = paste("optional named list with up to three elements,",
                                 "'coef', 'theta' and 'knots', specifying lower bounds",
                                 "for coefficients to be estimated.",
                                 "These must be finite and will be recycled if necessary to match",
                                 "`length(coefficients)`.")),
    defineParameter("nlminb.control", "numeric",
                    default = list(iter.max = 5e3L, eval.max = 5e3L),
                    desc = paste("optional list of control parameters to be passed to",
                                 "the `nlminb` optimizer. See `?nlminb`.")),
    defineParameter("plot", "logical", default = TRUE,
                    desc = "logical. Plot model fit against the data. Prediction interval"),
    defineParameter("start", "numeric, list", default = NULL,
                    desc = paste("optional starting values for the parameters to be estimated.
                                 Those are passed to `nlminb` and can be a single vector, or a list of vectors.",
                                 "In the latter case, only the best solution, that is,",
                                 "the one which minimizes the most the objective function, is kept.")),
    defineParameter("trace", "numeric", default = 1,
                    desc = paste("non-negative integer. If > 0, tracing information on the progress",
                                 "of the optimization are printed every `trace` iteration.",
                                 "If parallel computing is enable, nlminb trace logs are written",
                                 "into the working directory.",
                                 "Log files are prefixed with 'fireSense_IgnitionFit_trace'",
                                 "followed by the nodename (see ?Sys.info) and the subprocess pid.",
                                 "Default is 0, which turns off tracing.")),
    defineParameter("ub", "numeric", default = NULL,
                    desc = paste("optional named list with up to three elements,",
                                 "'coef', 'theta' and 'knots', specifying upper bounds",
                                 "for coefficients to be estimated.",
                                 "These must be finite and will be recycled if necessary to match",
                                 "`length(coefficients)`.")),
    defineParameter(".plots", "character", default = "screen",
                    desc = "See ?Plots. There are a few plots that are made within this module, if set."),
    defineParameter(".plotInitialTime", "numeric", default = start(sim),
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
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules,",
                                 "where stochasticity and time are not relevant."))
  ),
  inputObjects = expectsInput(
    objectName = "fireSense_ignitionCovariates",
    objectClass = "data.frame",
    desc = "One or more objects of class data.frame in which to look for variables present in the model formula.",
    sourceURL = NA_character_
  ),
  outputObjects = createsOutput(
    objectName = "fireSense_IgnitionFitted",
    objectClass = "fireSense_IgnitionFit",
    desc = "A fitted model object of class fireSense_IgnitionFit."
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_IgnitionFit = function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- current(sim)$moduleName

  switch(
    eventType,
    init = {
      sim <- frequencyFitInit(sim)

      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "run")

      if (!is.na(P(sim)$.saveInitialTime))
        sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, moduleName, "save", .last())
    },
    run = {
      sim <- frequencyFitRun(sim)

      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "run")
    },
    save = {
      sim <- frequencyFitSave(sim)

      if (!is.na(P(sim)$.saveInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, moduleName, "save", .last())
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )

  invisible(sim)
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
frequencyFitInit <- function(sim) {
  # Checking parameters
  stopifnot(P(sim)$trace >= 0)
  stopifnot(P(sim)$cores >= 1)
  stopifnot(P(sim)$iterDEoptim >= 1)
  stopifnot(P(sim)$iterNlminb >= 1)

  return(invisible(sim))
}

frequencyFitRun <- function(sim) {

  moduleName <- current(sim)$moduleName

  rescaleMDCFactor <- 1000
  fireSense_ignitionCovariates <- sim$fireSense_ignitionCovariates
  setDT(fireSense_ignitionCovariates)
  fireSense_ignitionCovariates[, MDC := MDC/rescaleMDCFactor]
  setDF(fireSense_ignitionCovariates)
  lb <- P(sim)$lb
  ub <- P(sim)$ub
  lb$knots <- lb$knots/rescaleMDCFactor
  ub$knots <- ub$knots/rescaleMDCFactor

  # Redo parameter bounds after rescale

  if (is.empty.model(as.formula(P(sim)$fireSense_ignitionFormula)))
    stop(moduleName, "> The formula describes an empty model.")

  # Remove rows of data with no cover and no ignitions
  whRowsHaveNoCover <- apply(as.data.frame(fireSense_ignitionCovariates)[,c(2,4:8)], 1, sum) == 0
  fireSense_ignitionCovariates <- fireSense_ignitionCovariates[!whRowsHaveNoCover,]

  # sim$fireSense_ignitionFormula <- paste0("ignitions ~ ",
  #                                         # "youngAge:MDC + ",
  #                                         "nonForest_highFlam:MDC + ",
  #                                         "nonForest_lowFlam:MDC + class2:MDC + class3:MDC + ",
  #                                         "youngAge:pw(MDC, k_YA) + nonForest_lowFlam:pw(MDC, k_NFLF) + ",
  #                                         # "nonForest_highFlam:pw(MDC, k_NFHF) + class2:pw(MDC, k_class2) + ",
  #                                         "class3:pw(MDC, k_class3) - 1")
  # params(sim)[[currentModule(sim)]]$fireSense_ignitionFormula <- sim$fireSense_ignitionFormula
  #
  fireSense_ignitionFormula <- as.formula(P(sim)$fireSense_ignitionFormula)
  terms <- terms.formula(fireSense_ignitionFormula, specials = "pw")

  if (attr(terms, "response")) {
    y <- fireSense_ignitionFormula[[2L]]
  } else {
    stop(moduleName, "> Incomplete formula, the LHS is missing.")
  }
  nx <- length(labels(terms)) + attr(terms, "intercept") ## Number of variables (covariates)
  allxy <- all.vars(terms)

  # Check the presence of at least one piecewise term
  hvPW <- !is.null(attr(terms, "specials")$pw)

  kLB <- kUB <- NULL

  if (hvPW) {
    objfun <- fireSenseUtils::.objFunIgnitionPW

    specialsInd <- which(unlist(lapply(attr(terms,"variables"), is.call)))
    specialsCalls <- attr(terms,"variables")[specialsInd]

    ## Extract the names of the knots (breakpoints)
    ## Alternative way: all.vars(terms)[!all.vars(terms) %in% rownames(attr(terms,"factors"))]
    specialsTerms <- lapply(specialsCalls, function(specialsCall) {
      if (specialsCall[[1L]] == "pw") {
        specialsCall[[1L]] <- quote(extractSpecial)
        eval(specialsCall)
      }
    }
    )

    specialsTerms <- specialsTerms[!unlist(lapply(specialsTerms, is.null))]

    kNames <- sapply(specialsTerms, "[[", "knot")

    if (anyDuplicated(kNames))
      stop(moduleName, "> Knot's names are not unique.")

    nk <- length(kNames)
    allx <- allxy[!allxy %in% c(y, kNames)]

    missing <- !allxy[!allxy %in% kNames] %in% ls(fireSense_ignitionCovariates, all.names = TRUE)

    if (s <- sum(missing))
      stop(moduleName, "> '", allxy[!allxy %in% kNames][missing][1L], "'",
           if (s > 1) paste0(" (and ", s - 1L, " other", if (s > 2) "s", ")"),
           " not found in data objects nor in the simList environment.")

    ## Covariates that have a breakpoint
    pwVarNames <- sapply(specialsTerms, "[[", "variable", simplify = FALSE)

    kUB <- if (is.null(ub$knots)) {
      lapply(pwVarNames, function(x) max(if (is(x, "AsIs")) x else fireSense_ignitionCovariates[[x]])) %>% unlist()
    } else {
      rep_len(ub$knots, nk) ## User-defined bounds (recycled if necessary)
    }

    kLB <- if (is.null(lb$knots)) {
      lapply(pwVarNames, function(x) min(if (is(x, "AsIs")) x else fireSense_ignitionCovariates[[x]])) %>% unlist()
    } else {
      rep_len(lb$knots, nk) ## User-defined bounds (recycled if necessary)
    }

    knots <- mapply(
      kNames,
      z = pwVarNames,
      FUN = function(w, z)
        #fireSense_ignitionCovariates[[w]] <-
        mean(if (is(z, "AsIs")) z else fireSense_ignitionCovariates[[z]]),
      SIMPLIFY = FALSE
    )
    fireSense_ignitionCovariates <- data.frame(fireSense_ignitionCovariates, knots)

    updateKnotExpr <- parse(text = paste0("mod_env[[\"", kNames, "\"]] = params[", (nx + 1L):(nx + nk), "]", collapse = "; "))
  } else {
    missing <- !allxy %in% ls(fireSense_ignitionCovariates, all.names = TRUE)

    if (s <- sum(missing))
      stop(moduleName, "> '", allxy[missing][1L], "'",
           if (s > 1) paste0(" (and ", s - 1L, " other", if (s > 2) "s", ")"),
           " not found in data objects nor in the simList environment.")

    allx <- allxy[allxy != y]
    objfun <- fireSenseUtils::.objFunIgnition
    nk <- 0L
  }

  family <- P(sim)$family

  if (is.language(family)) family <- eval(family)

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
    family <- tryCatch(
      family(),
      error = function(e) family(
        theta = suppressWarnings(
          glm.nb(formula = fireSense_ignitionFormula,
                 y = FALSE,
                 model = FALSE,
                 data = fireSense_ignitionCovariates)[["theta"]]
        )
      )
    )
  } else if (is.function(family)) {
    family <- family()
  }

  # Check ff family is negative.binomial
  isFamilyNB <- grepl(x = family$family, pattern = "Negative Binomial")

  # Extract starting value for theta
  if (isFamilyNB) {
    theta <- as.numeric(gsub(x = regmatches(family$family, gregexpr("\\(.*?\\)", family$family))[[1L]],
                             pattern = "[\\(\\)]", replacement = ""))
  }

  linkinv <- family$linkinv

  mm <- model.matrix(object = terms, data = fireSense_ignitionCovariates)

  # Does the model formula contain an offset?
  model_offset <- model.offset(model.frame(fireSense_ignitionFormula, fireSense_ignitionCovariates))
  offset <- if (is.null(model_offset)) 0 else model_offset

  ## Define the scaling matrix. This is used later in the optimization process
  ## to rescale parameter values between 0 and 1, i.e. put all variables on the same scale.
  n <- nx + nk
  if (isFamilyNB) n <- n + 1L

  sm <- matrix(0, n, n)
  diag(sm) <- 1

  ## Define parameter bounds automatically if they are not supplied by user
  ## First defined the bounds for DEoptim, the first optimizer

  ## Upper bounds
  DEoptimUB <- c(
    if (is.null(ub[["coef"]])) {
      ## Automatically estimate an upper boundary for each parameter
      (suppressWarnings(
        tryCatch(
          glm(
            formula = fireSense_ignitionFormula,
            y = FALSE,
            model = FALSE,
            data = fireSense_ignitionCovariates,
            family = poisson(link = family$link)
          ),
          error = function(e) stop(
            moduleName, "> Automated estimation of upper bounds",
            " (coefs) failed, please set the 'coef' element of ",
            "the 'ub' parameter."
          )
        )
      ) %>% coef() %>% oom(.)) * 10L -> ub

      if (anyNA(ub)) {
        stop(moduleName, "> Automated estimation of upper bounds (coefs) failed, ",
             "please set the 'coef' element of the 'ub' parameter.")
      } else {
        ub
      }
    } else {
      rep_len(ub[["coef"]], nx) ## User-defined bounds (recycled if necessary)
    },
    kUB
  )

  ## Lower bounds
  DEoptimLB <- c({
    switch(family$link,
           log = {
             if (is.null(lb[["coef"]])) {
               -DEoptimUB[1L:nx] ## Automatically estimate a lower boundary for each parameter
             } else {
               rep_len(lb[["coef"]], nx) ## User-defined bounds (recycled if necessary)
             }

           }, identity = {
             if (is.null(lb[["coef"]])) {
               rep_len(1e-16, nx) ## Ensure non-negativity
             } else {
               rep_len(lb[["coef"]], nx) ## User-defined bounds (recycled if necessary)
             }
           }, stop(moduleName, "> Link function ", family$link, " is not supported."))
  }, kLB)

  ## If negative.binomial family needs to add bounds for theta parameter
  if (isFamilyNB) {
    DEoptimUB <- c(DEoptimUB, if (is.null(ub$t)) 2L * theta else ub$t)
    DEoptimLB <- c(DEoptimLB, if (is.null(lb$t)) 1e-16 else lb$t) ## Enfore non-negativity
  }

  ## Then, define lower and upper bounds for the second optimizer (nlminb)
  ## Upper bounds
  nlminbUB <- DEoptimUB
  if (is.null(ub[["coef"]]))
    nlminbUB[1:nx] <- rep_len(Inf, nx)

  ## Lower bounds
  nlminbLB <- if (is.null(lb[["coef"]])) {
    c(switch(family$link,
             log = rep_len(-Inf, nx),       ## log-link, default: -Inf for terms and 0 for breakpoints/knots
             identity = rep_len(1e-16, nx)) ## identity link, default: enforce non-negativity
      , kLB)
  } else {
    DEoptimLB ## User-defined lower bounds for parameters to be estimated
  }

  ## If negative.binomial family add bounds for the theta parameter
  if (isFamilyNB && is.null(lb$t)) {
    nlminbLB <- c(nlminbLB, 1e-16) ## Enforce non-negativity
  } else if (isFamilyNB) {
    nlminbLB <- c(nlminbLB, lb$t)
  }

  ## Define the log-likelihood function (objective function)
  nll <- switch(family$family,
                poisson = parse(text = paste0("-sum(dpois(x=", y, ", lambda = mu, log = TRUE))")),
                parse(text = paste0("-sum(dnbinom(x=", y, ", mu = mu, size = params[length(params)], log = TRUE))")))

  trace <- P(sim)$trace
  if (P(sim)$cores > 1) {
    if (.Platform$OS.type == "unix") {
      mkCluster <- parallel::makeForkCluster
    } else {
      #mkCluster <- parallel::makePSOCKcluster ## TODO: this attaches `snow` and breaks the module
      ## see warning: https://www.rdocumentation.org/packages/secr/versions/4.3.3/topics/Parallel
    }

    message("Creating cluster")
    cl <- mkCluster(P(sim)$cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library("MASS"))
    # assign("mod_env", fireSense_ignitionCovariates, envir = .GlobalEnv)
    # clusterExport(cl, varlist = list("mod_env"), envir = environment()) # it is faster to "get" it internally
  }

  ## If starting values are not supplied
  if (is.null(P(sim)$start)) {
    ## First optimizer, get rough estimates of the parameter values
    ## Use these estimates to compute the order of magnitude of these parameters

    control <- list(itermax = P(sim)$iterDEoptim, trace = P(sim)$trace)
    if (P(sim)$cores > 1) {
      control$cluster <- cl
    }

    # control$cluster <- NULL #when debugging DEoptim
    if (hvPW) {
      DEout <- Cache(DEoptim, objfun, lower = DEoptimLB, upper = DEoptimUB,
                     control = do.call("DEoptim.control", control),
                     formula = P(sim)$fireSense_ignitionFormula,
                     mod_env = fireSense_ignitionCovariates,
                     linkinv = linkinv, nll = nll, sm = sm, nx = nx,
                     offset = offset, updateKnotExpr = updateKnotExpr,
                     userTags = c("ignitionFit", "DEoptim"))
      DEoptimBestMem <- DEout %>%  `[[`("optim") %>% `[[`("bestmem")
    } else {
      DEoptimBestMem <- Cache(DEoptim, objfun, lower = DEoptimLB, upper = DEoptimUB,
                              control = do.call("DEoptim.control", control),
                              formula = P(sim)$fireSense_ignitionFormula,
                              mod_env = fireSense_ignitionCovariates,
                              linkinv = linkinv, nll = nll, sm = sm, nx = nx, mm = mm,
                              offset = offset,
                              userTags = c("ignitionFit", "DEoptim")) %>%
        `[[`("optim") %>% `[[`("bestmem")
    }

    if (FALSE) { # THIS IS ELIOT BLOCKING THIS -- I DON"T THINK IT IS NECESSARY
      ## Update scaling matrix
      diag(sm) <- oom(DEoptimBestMem)

      #TODO: I believe this is a mistake, and the length(nlminUB) subset should be replaced with
      ## Update of the lower and upper bounds of the coefficients based on the scaling matrix
      nlminbLB <- nlminbLB / diag(sm)
      nlminbUB <- nlminbUB / diag(sm)

      #nlminbLB[c(1:nx, length(nlminbLB))] <- nlminbLB[c(1:nx, length(nlminbLB))] / diag(sm)[c(1:nx, length(nlminbLB))]
      #nlminbUB[c(1:nx, length(nlminbUB))] <- nlminbUB[c(1:nx, length(nlminbUB))] / diag(sm)[c(1:nx, length(nlminbUB))]

      ## Update of the lower and upper bounds for the knots based on the scaling matrix
      #TODO: fix the partial matching

      if (hvPW) {
        kLB <- DEoptimLB[(nx + 1L):(nx + nk)] / diag(sm)[(nx + 1L):(nx + nk)]
        nlminbLB[(nx + 1L):(nx + nk)] <- if (is.null(lb$knots)) kLB else pmax(kLB, lb$knots)

        kUB <- DEoptimUB[(nx + 1L):(nx + nk)] / diag(sm)[(nx + 1L):(nx + nk)]
        nlminbUB[(nx + 1L):(nx + nk)] <- if (is.null(ub$knots)) kUB else pmin(kUB, ub$knots)
      }

      getRandomStarts <- function(.) {
       pmin(pmax(rnorm(length(DEoptimBestMem), 0L, 2L) / 40 +
                   unname(DEoptimBestMem / oom(DEoptimBestMem)), nlminbLB), nlminbUB)
      }
      start <- c(lapply(1:P(sim)$iterNlminb, getRandomStarts),
                list(unname(DEoptimBestMem / oom(DEoptimBestMem))))
    }

    # This is ELIOT's change March 5, 2021 -- simpler -- use more Information from DEoptim
    start <- split(DEout$member$pop, seq(NROW(DEout$member$pop)))

    if (TRUE) { # This allows a faster estimate, but fewer different starting values
      pop <- if(exists("cl", inherits = FALSE)) length(cl) else 1
      message("Subsetting for nlminb -- only running ", pop, " of the DEoptim pops")
      subset <- sample(NROW(start), size = pop)
      start <- start[subset]
    }
  } else {
    start <- if (is.list(P(sim)$start)) {
      diag(sm) <- lapply(P(sim)$start, oom) %>%
        do.call("rbind", .) %>%
        apply(2, function(x) as.numeric(names(base::which.max(table(x)))))

      lapply(P(sim)$start, function(x) x / diag(sm))
    } else {
      diag(sm) <- oom(P(sim)$start)
      P(sim)$start / diag(sm)
    }
  }

  outNlminb <- if (is.list(start)) {

    if (P(sim)$cores > 1) {
      outputPath <- outputPath(sim)
      basePattern <- paste(moduleName, Sys.info()[["nodename"]], format(Sys.time(), "%Y%m%d"), "trace", sep = "_")

      if (trace) {
        clusterExport(cl, c("outputPath", "basePattern"), envir = environment())

        clusterEvalQ(
          cl,
          sink(file.path(outputPath, paste0(basePattern, ".", Sys.getpid())))
        )
        pids <- unlist(clusterEvalQ(cl, Sys.getpid()))
        message("These are the pids of the spawned cores: ")
        dput(pids)

      }
      message("Starting nlminb ... ")
      out <- Cache(clusterApplyLB, cl = cl, x = start, fun = objNlminb, objective = objfun,
                   lower = nlminbLB, upper = nlminbUB, hvPW = hvPW,
                   linkinv = linkinv, nll = nll, sm = sm, nx = nx, mm = mm, #TODO mm may not be required with PW...
                   mod_env = fireSense_ignitionCovariates, offset = offset,
                   formula = P(sim)$fireSense_ignitionFormula,
                   omitArgs = c("x"), # don't need to know the random sample... the mm is enough
                   updateKnotExpr = updateKnotExpr, # cacheId = "e016b5d728ed2b6a",
                   control = c(P(sim)$nlminb.control, list(trace = trace)))

      if (FALSE) { # THIS SECTION ALLOWS MANUAL READING OF LOG FILES
        #  MUST MANUALLY IDENTIFY THE PIDS
        if (exists("pids")) {
          aa <- lapply(paste0("~/GitHub/WBI_fireSense/outputs/AB/fireSense_IgnitionFit_spades189_20210308_trace.", pids), readLines)
          bb <- lapply(aa, function(bb) {
            bb <- gsub("^ +", "", bb)
            vals <- strsplit(bb, split = ":* +")
            cc <- do.call(rbind, lapply(vals, as.numeric))
            wh <- unique(c(which(cc[, 1] == 0) - 1, NROW(cc)))
            wh <- setdiff(wh, 0)
            cc[wh,, drop = FALSE]
          })
          cc <- do.call(rbind, bb)
          dd <- head(data.table::as.data.table(cc[order(cc[, 2]),]), 20)
          colnms <- c("IterationNumStopped", "ObjFunValue", attr(terms, "term.labels"),
                      unlist(lapply(updateKnotExpr, function(x) x[[2]][[3]])), "NB_theta")
          (data.table::setnames(dd, colnms))
        }
      }
      message("... Done")

      if (trace) clusterEvalQ(cl, sink())
    } else {
      warning("This is not tested by Eliot as of March 4, 2021; please set parameter: cores > 1")
      out <- Cache(lapply, start[1], objNlminb, objective = objfun, lower = nlminbLB, upper = nlminbUB, hvPW = hvPW,
                    linkinv = linkinv, nll = nll, sm = sm, nx = nx, mm = mm, #TODO mm may not be required with PW...
                    mod_env = fireSense_ignitionCovariates, offset = offset,
                    formula = P(sim)$fireSense_ignitionFormula,
                    updateKnotExpr = updateKnotExpr, omitArgs = "X",
                    control = c(P(sim)$nlminb.control, list(trace = min(6, trace * 3))))
    }

    out
  } else {
    warning("This is not tested by Eliot as of March 4, 2021; please run more than one start")
    list(objNlminb(start, objfun, nlminbLB, nlminbUB, c(P(sim)$nlminb.control, list(trace = trace))))
  }

  ## Select best minimum amongst all trials
  outBest <- outNlminb[[which.min(sapply(outNlminb, "[[", "objective"))]]

  ## Compute the standard errors around the estimates
  message("Computing hessian by numerical derivation")
  if (hvPW) {
    hess <- Cache(numDeriv::hessian, func = objfun, x = outBest$par,
                  mod_env = fireSense_ignitionCovariates,
                  linkinv = linkinv, nll = nll,
                  sm = sm, nx = nx, updateKnotExpr = updateKnotExpr,
                  formula = P(sim)$fireSense_ignitionFormula,
                  offset = offset,
                  userTags = c("fireSense_ignitionFit", "hessian"))
  } else {
    hess <- Cache(numDeriv::hessian, func = objfun, x = outBest$par,
                  mod_env = fireSense_ignitionCovariates,
                  linkinv = linkinv, nll = nll,
                  sm = sm, nx = nx, mm = mm,
                  offset = offset,
                  userTags = c("fireSense_ignitionFit", "hessian"))
  }

  solvedHess <- tryCatch(solve(hess), error = function(e) NA)
  se <- suppressWarnings(tryCatch(drop(sqrt(diag(solvedHess)) %*% sm), error = function(e) NA))
  message("... Done")

  if (!exists("outBest", inherits = FALSE)) {
    best <- drop(as.matrix(dd[1, -(1:2)]))
  } else {
    best <- outBest$par
    colnms <- c(attr(terms, "term.labels"),
                unlist(lapply(updateKnotExpr, function(x) x[[2]][[3]])), "NB_theta")
    names(best) <- colnms
  }

  if (anyPlotting(P(sim)$.plots)) {

    ndLong <- pwPlotData(bestParams = best,
                         ses = se, solvedHess = solvedHess,
                         formula = P(sim)$fireSense_ignitionFormula,
                         xColName = "MDC", nx = nx, offset = offset, linkinv = linkinv)

    Plots(data = ndLong, fn = pwPlot,
          ggTitle =  paste0(basename(outputPath(sim)), " fireSense IgnitionFit"),
          filename = "IgnitionRatePer100km2")#, types = "screen", .plotInitialTime = time(sim))
   #TODO unresolved bug in Plot triggered by spaces
  }


  convergence <- TRUE

  if (outBest$convergence || anyNA(se)) {
    tooClose <- 0.00001
    closeToBounds <- abs(drop((best - DEoptimLB)/ (DEoptimUB - DEoptimLB))) < tooClose |
      best < tooClose
    ctb <- data.table(term = names(closeToBounds), best = best,
                      upperBoundary = DEoptimUB,
                      lowerBoundary = DEoptimLB,
                      closeToBounds = closeToBounds)[closeToBounds]
    possTerms <- attr(terms, "term.labels")[1:nx][!closeToBounds[1:nx]]
    possForm <- paste0(terms[[2]], " ~ ", paste(possTerms, collapse = " + "), " -1")


    message("--------------------------------------------------")
    message("It is possible that parameters are too close to their boundary values (or zero). ",
            "The following are within ",tooClose*100,"% of their boundary and removing them ",
            "from sim$fireSense_ignitionFormula may help with convergence or invertability... e.g.")
    message("sim$fireSense_ignitionFormula <- \"", possForm, "\"")
    messageDF(ctb)
    message("It may also help to use Ben Bolker's approximation: sqrt(1/diag(hess)) mentioned here:")
    message("https://cran.r-project.org/web/packages/bbmle/vignettes/mle2.pdf")
    message("If there are Inf values, that indicates variables to remove as they have",
            "infinite variance at the solution")

    if (!isFALSE(P(sim)$autoRefit))  {
      outRL <- if (isTRUE(P(sim)$autoRefit)) {
        message("Automatically refitting with simpler model becaues P(sim)$autoRefit is TRUE")
        "y"
      } else if (isFALSE(P(sim)$autoRefit)) {
        "n"
      } else {
        readline("Would you like to restart this IgnitionFit event with that new formula (Y or N or interactive)? ")
      }
      if (identical(tolower(outRL), "y")) {
        params(sim)[[currentModule(sim)]]$fireSense_ignitionFormula <- possForm
        sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "run", eventPriority = 1)
      } else if (isTRUE(startsWith(outRL, "i"))) {
        browser()
      }
    }

    if (outBest$convergence) {
      convergence <- FALSE
      convergDiagnostic <- paste0("nlminb optimizer did not converge (", outBest$message, ")")
      warning(moduleName, "> ", convergDiagnostic, immediate. = TRUE)
    } else if (anyNA(se)) {
      ## Negative values in the Hessian matrix suggest that the algorithm did not converge

      convergence <- FALSE
      convergDiagnostic <- "nlminb optimizer reached relative convergence, saddle point?"
      # warning(moduleName, "> ", convergDiagnostic, immediate. = TRUE)
    }
  } else {
    convergDiagnostic <- outBest$message
  }

  ## Parameters scaling: Revert back estimated coefficients to their original scale
  outBest$par <- drop(outBest$par %*% sm)

  l <- list(formula = as.formula(fireSense_ignitionFormula),
            family = family,
            coef = setNames(outBest$par[1:nx], colnames(mm)),
            coef.se = setNames(se[1:nx], colnames(mm)),
            LL = -outBest$objective,
            AIC = 2 * length(outBest$par) + 2 * outBest$objective,
            convergence = convergence,
            convergenceDiagnostic = convergDiagnostic,
            rescales = list(MDC = paste0("MDC / ", rescaleMDCFactor)))

  if (hvPW) {
    l$knots <- setNames(outBest$par[(nx + 1L):(nx + nk)], kNames)
    l$knots.se <- setNames(se[(nx + 1L):(nx + nk)], kNames)
  }

  if (isFamilyNB) {
    ## Update the NB family template with the estimated theta
    l$family <- MASS::negative.binomial(theta = outBest$par[length(outBest$par)], link = family$link)
    l$theta <- outBest$par[length(outBest$par)]
    l$theta.se <- se[length(se)]
  }

  sim$fireSense_IgnitionFitted <- l
  class(sim$fireSense_IgnitionFitted) <- "fireSense_IgnitionFit"

  invisible(sim)
}

frequencyFitSave <- function(sim) {
  timeUnit <- timeunit(sim)
  currentTime <- time(sim, timeUnit)

  saveRDS(
    sim$fireSense_IgnitionFitted,
    file = file.path(paths(sim)$out, paste0("fireSense_IgnitionFitted_", timeUnit, currentTime, ".rds"))
  )

  invisible(sim)
}

pwPlot <- function(d, ggTitle)  {
  gg <- ggplot(d,  aes(x=MDC, y=mu, group=Type, color=Type)) +
    geom_line()
  if (!anyNA(d$lci))
    gg <- gg + geom_smooth(aes(ymin = lci, ymax = uci), stat = "identity")
  gg <- gg +
    labs(y = "Igntion rate per 100 km2", title = ggTitle) +
    theme_bw()
}

pwPlotData <- function(bestParams, formula, xColName = "MDC", nx, offset, linkinv,
                       solvedHess, ses) {

  newDat <- expand.grid(MDC = 1:250/1000, youngAge = 0:1, nonForest_highFlam = 0:1,
                        nonForest_lowFlam = 0:1,
                        class2 = 0:1, class3 = 0:1)
  cns <- setdiff(colnames(newDat), xColName)

  # bestParams <- drop(as.matrix(dd[1, -(1:2)]))

  keep <- grep("^k_", names(bestParams), value = TRUE)
  newDat <- data.table(newDat, t(bestParams[keep]))
  for (cn in cns) {
    set(newDat, which(newDat[[cn]] == 1), setdiff(cns, cn), 0)
  }
  newDat <- unique(newDat)
  keepers <- apply(newDat[, ..cns], 1, function(x) sum(x) > 0)
  newDat <- newDat[keepers]

  mm <- model.matrix(as.formula(formula)[-2], newDat)
  mu <- drop(mm %*% bestParams[1:nx]) + offset
  if (!all(is.na(solvedHess))) {
    # https://biologyforfun.wordpress.com/2015/06/17/confidence-intervals-for-prediction-in-glmms/
    pvar1 <- diag(mm %*% tcrossprod(solvedHess[1:nx, 1:nx], mm))

    uci <- mu + sqrt(pvar1) * 1.96 + offset
    lci <- mu - sqrt(pvar1) * 1.96 + offset
  } else {
    uci <- NA_real_
    lci <- NA_real_
  }

  ## link implementation
  mu <- linkinv(mu)
  newDat$mu <- mu
  newDat <- newDat[, !..keep]
  newDat[, MDC := MDC * 1000]
  newDat[, `:=`(lci  = lci, uci = uci)]
  setkeyv(newDat, "MDC")
  ndLong <- data.table::melt(newDat, measure.vars = 2:6, variable.name = "Type")
  ndLong <- ndLong[value > 0]
}
