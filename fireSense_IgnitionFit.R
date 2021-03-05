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
  reqdPkgs = list("DEoptim", "dplyr", "MASS", "magrittr", "numDeriv", "parallel", "pemisc",
                  "PredictiveEcology/fireSenseUtils@development (>=0.0.4.9044)"),
  parameters = rbind(
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
    defineParameter("iterDEoptim", "integer", default = 2000,
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

  setDT(sim$fireSense_ignitionCovariates)
  sim$fireSense_ignitionCovariates[, MDC := MDC/1000]
  setDF(sim$fireSense_ignitionCovariates)
  # Redo parameter bounds after rescale
  params(sim)[[currentModule(sim)]]$lb$knots <- quantile(sim$fireSense_ignitionCovariates$MDC, probs = 0.1)
  params(sim)[[currentModule(sim)]]$ub$knots <- quantile(sim$fireSense_ignitionCovariates$MDC, probs = 0.9)

  if (is.empty.model(as.formula(P(sim)$fireSense_ignitionFormula)))
    stop(moduleName, "> The formula describes an empty model.")

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

    missing <- !allxy[!allxy %in% kNames] %in% ls(sim$fireSense_ignitionCovariates, all.names = TRUE)

    if (s <- sum(missing))
      stop(moduleName, "> '", allxy[!allxy %in% kNames][missing][1L], "'",
           if (s > 1) paste0(" (and ", s - 1L, " other", if (s > 2) "s", ")"),
           " not found in data objects nor in the simList environment.")

    ## Covariates that have a breakpoint
    pwVarNames <- sapply(specialsTerms, "[[", "variable", simplify = FALSE)

    kUB <- if (is.null(P(sim)$ub$k)) {
      lapply(pwVarNames, function(x) max(if (is(x, "AsIs")) x else sim$fireSense_ignitionCovariates[[x]])) %>% unlist()
    } else {
      rep_len(P(sim)$ub$k, nk) ## User-defined bounds (recycled if necessary)
    }

    kLB <- if (is.null(P(sim)$lb$k)) {
      lapply(pwVarNames, function(x) min(if (is(x, "AsIs")) x else sim$fireSense_ignitionCovariates[[x]])) %>% unlist()
    } else {
      rep_len(P(sim)$lb$k, nk) ## User-defined bounds (recycled if necessary)
    }

    invisible(
      mapply(
        kNames,
        z = pwVarNames,
        FUN = function(w, z)
          sim$fireSense_ignitionCovariates[[w]] <- mean(if (is(z, "AsIs")) z else sim$fireSense_ignitionCovariates[[z]]),
        SIMPLIFY = FALSE
      )
    )

    updateKnotExpr <- parse(text = paste0("mod_env[[\"", kNames, "\"]] = params[", (nx + 1L):(nx + nk), "]", collapse = "; "))
  } else {
    missing <- !allxy %in% ls(sim$fireSense_ignitionCovariates, all.names = TRUE)

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
                 data = sim$fireSense_ignitionCovariates)[["theta"]]
        )
      )
    )
  } else if (is.function(family)) {
    family <- family()
  }

  browser()
  # Check ff family is negative.binomial
  isFamilyNB <- grepl(x = family$family, pattern = "Negative Binomial")

  # Extract starting value for theta
  if (isFamilyNB) {
    theta <- as.numeric(gsub(x = regmatches(family$family, gregexpr("\\(.*?\\)", family$family))[[1L]],
                             pattern = "[\\(\\)]", replacement = ""))
  }

  linkinv <- family$linkinv

  mm <- model.matrix(object = terms, data = sim$fireSense_ignitionCovariates)

  # Does the model formula contain an offset?
  model_offset <- model.offset(model.frame(fireSense_ignitionFormula, sim$fireSense_ignitionCovariates))
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
    if (is.null(P(sim)$ub[["coef"]])) {
      ## Automatically estimate an upper boundary for each parameter
      (suppressWarnings(
        tryCatch(
          glm(
            formula = fireSense_ignitionFormula,
            y = FALSE,
            model = FALSE,
            data = sim$fireSense_ignitionCovariates,
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
      rep_len(P(sim)$ub[["coef"]], nx) ## User-defined bounds (recycled if necessary)
    },
    kUB
  )

  ## Lower bounds
  DEoptimLB <- c({
    switch(family$link,
           log = {
             if (is.null(P(sim)$lb[["coef"]])) {
               -DEoptimUB[1L:nx] ## Automatically estimate a lower boundary for each parameter
             } else {
               rep_len(P(sim)$lb[["coef"]], nx) ## User-defined bounds (recycled if necessary)
             }

           }, identity = {
             if (is.null(P(sim)$lb[["coef"]])) {
               rep_len(1e-16, nx) ## Ensure non-negativity
             } else {
               rep_len(P(sim)$lb[["coef"]], nx) ## User-defined bounds (recycled if necessary)
             }
           }, stop(moduleName, "> Link function ", family$link, " is not supported."))
  }, kLB)

  ## If negative.binomial family needs to add bounds for theta parameter
  if (isFamilyNB) {
    DEoptimUB <- c(DEoptimUB, if (is.null(P(sim)$ub$t)) 2L * theta else P(sim)$ub$t)
    DEoptimLB <- c(DEoptimLB, if (is.null(P(sim)$lb$t)) 1e-16 else P(sim)$lb$t) ## Enfore non-negativity
  }

  ## Then, define lower and upper bounds for the second optimizer (nlminb)
  ## Upper bounds
  nlminbUB <- DEoptimUB
  if (is.null(P(sim)$ub[["coef"]]))
    nlminbUB[1:nx] <- rep_len(Inf, nx)

  ## Lower bounds
  nlminbLB <- if (is.null(P(sim)$lb[["coef"]])) {
    c(switch(family$link,
             log = rep_len(-Inf, nx),       ## log-link, default: -Inf for terms and 0 for breakpoints/knots
             identity = rep_len(1e-16, nx)) ## identity link, default: enforce non-negativity
      , kLB)
  } else {
    DEoptimLB ## User-defined lower bounds for parameters to be estimated
  }

  ## If negative.binomial family add bounds for the theta parameter
  if (isFamilyNB && is.null(P(sim)$lb$t)) {
    nlminbLB <- c(nlminbLB, 1e-16) ## Enforce non-negativity
  } else if (isFamilyNB) {
    nlminbLB <- c(nlminbLB, P(sim)$lb$t)
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
      mkCluster <- parallel::makePSOCKcluster
    }

    message("Creating cluster")
    cl <- mkCluster(P(sim)$cores)
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, library("MASS"))
    # assign("mod_env", sim$fireSense_ignitionCovariates, envir = .GlobalEnv)
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

      control$itermax <- 300
      browser()
      DEout <- Cache(DEoptim, objfun, lower = DEoptimLB, upper = DEoptimUB,
                              control = do.call("DEoptim.control", control),
                              formula = P(sim)$fireSense_ignitionFormula,
                     mod_env = sim$fireSense_ignitionCovariates,
                     linkinv = linkinv, nll = nll, sm = sm, nx = nx,
                     offset = offset, updateKnotExpr = updateKnotExpr,
                     userTags = c("ignitionFit", "DEoptim"))
      DEoptimBestMem <- DEout %>%  `[[`("optim") %>% `[[`("bestmem")
    } else {
      DEoptimBestMem <- Cache(DEoptim, objfun, lower = DEoptimLB, upper = DEoptimUB,
                              control = do.call("DEoptim.control", control),
                              formula = P(sim)$fireSense_ignitionFormula,
                              mod_env = sim$fireSense_ignitionCovariates,
                              linkinv = linkinv, nll = nll, sm = sm, nx = nx, mm = mm,
                              offset = offset,
                              userTags = c("ignitionFit", "DEoptim")) %>%
        `[[`("optim") %>% `[[`("bestmem")
    }

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
      nlminbLB[(nx + 1L):(nx + nk)] <- if (is.null(P(sim)$lb$k)) kLB else pmax(kLB, P(sim)$lb$k)

      kUB <- DEoptimUB[(nx + 1L):(nx + nk)] / diag(sm)[(nx + 1L):(nx + nk)]
      nlminbUB[(nx + 1L):(nx + nk)] <- if (is.null(P(sim)$ub$k)) kUB else pmin(kUB, P(sim)$ub$k)
    }

    start <- split(DEout$member$pop, seq(NROW(DEout$member$pop)))
    #getRandomStarts <- function(.) {
    #  pmin(pmax(rnorm(length(DEoptimBestMem), 0L, 2L) / 40 +
    #              unname(DEoptimBestMem / oom(DEoptimBestMem)), nlminbLB), nlminbUB)
    #}
    #start <- c(lapply(1:P(sim)$iterNlminb, getRandomStarts),
    #           list(unname(DEoptimBestMem / oom(DEoptimBestMem))))
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

  out <- if (is.list(start)) {

    if (P(sim)$cores > 1) {
      outputPath <- outputPath(sim)
      basePattern <- paste(moduleName, Sys.info()[["nodename"]], format(Sys.time(), "%Y%m%d"), "trace", sep = "_")

      if (trace) {
        clusterExport(cl, c("outputPath", "basePattern"), envir = environment())

        clusterEvalQ(
          cl,
          sink(file.path(outputPath, paste0(basePattern, ".", Sys.getpid())))
        )
      }
      message("Starting nlminb ... ")
      out <- clusterApplyLB(cl = cl, x = start, fun = objNlminb, objective = objfun,
                            lower = nlminbLB, upper = nlminbUB, hvPW = hvPW,
                            linkinv = linkinv, nll = nll, sm = sm, nx = nx, mm = mm, #TODO mm may not be required with PW...
                            mod_env = sim$fireSense_ignitionCovariates, offset = offset,
                            formula = P(sim)$fireSense_ignitionFormula,
                            updateKnotExpr = updateKnotExpr,
                            control = c(P(sim)$nlminb.control, list(trace = trace)))
      if (FALSE) {
        pids <- 1784000 + 307:320; aa <- lapply(paste0("~/GitHub/WBI_fireSense/outputs/NT/fireSense_IgnitionFit_spades189_20210304_trace.", pids), readLines)
        bb <- lapply(aa, function(bb) {
          cc <- do.call(rbind, lapply(strsplit(bb, split = ":* +"), as.numeric))
          wh <- which(cc[, 2] == 0) - 1
          wh <- setdiff(wh, 0)
          cc[wh,, drop = FALSE]
        })
        cc <- do.call(rbind, bb)
        dd <- as.data.table(cc[order(cc[, 2]),])
      }
      message("... Done")

      if (trace) clusterEvalQ(cl, sink())
    } else {
      DEoptimBestMem <- c(0.035547,    0.005960,    0.001118,    0.318472,    0.227213,    0.803050,    0.106184,    0.019238,    0.074035,    0.306645,    0.112640,    0.120402,    0.142217,    0.138744,    0.109269)

      # DEoptimBestMem <- c(0.000010,    0.000007,    0.000003,    0.000289,    0.000284,    0.001340,    0.000091,    0.000020,    0.000034,    0.000269,  123.202095,  124.944181,  172.975829,  172.509741,  138.399639)
      # DEoptimBestMem <- c(0.000007,    0.000009,    0.000001,    0.000292,    0.000282,    0.001348,    0.000094,    0.000029,    0.000023,    0.000260,  124.729139,  122.971900,  172.988124,  172.945087,  135.848640)
      start <- list(DEoptimBestMem)
      out <- lapply(start, objNlminb, objective = objfun, lower = nlminbLB, upper = nlminbUB, hvPW = hvPW,
                    linkinv = linkinv, nll = nll, sm = sm, nx = nx, mm = mm, #TODO mm may not be required with PW...
                    mod_env = sim$fireSense_ignitionCovariates, offset = offset,
                    formula = P(sim)$fireSense_ignitionFormula,
                    updateKnotExpr = updateKnotExpr,
                    control = c(P(sim)$nlminb.control, list(trace = trace)))
    }

    ## Select best minimum amongst all trials
    out[[which.min(sapply(out, "[[", "objective"))]]
  } else {
    objNlminb(start, objfun, nlminbLB, nlminbUB, c(P(sim)$nlminb.control, list(trace = trace)))
  }

  ## Compute the standard errors around the estimates
  if (hvPW) {
    hess <- Cache(numDeriv::hessian, func = objfun, x = out$par,
                  mod_env = sim$fireSense_ignitionCovariates,
                  linkinv = linkinv, nll = nll,
                  sm = sm, nx = nx, updateKnotExpr = updateKnotExpr,
                  formula = P(sim)$fireSense_ignitionFormula,
                  offset = offset,
                  userTags = c("fireSense_ignitionFit", "hessian"))
  } else {
    hess <- Cache(numDeriv::hessian, func = objfun, x = out$par,
                  mod_env = sim$fireSense_ignitionCovariates,
                  linkinv = linkinv, nll = nll,
                  sm = sm, nx = nx, mm = mm,
                  offset = offset,
                  userTags = c("fireSense_ignitionFit", "hessian"))
  }

  se <- suppressWarnings(tryCatch(drop(sqrt(diag(solve(hess))) %*% sm), error = function(e) NA))

  convergence <- TRUE

  if (out$convergence) {
    convergence <- FALSE
    convergDiagnostic <- paste0("nlminb optimizer did not converge (", out$message, ")")
    warning(moduleName, "> ", convergDiagnostic, immediate. = TRUE)
  } else if (anyNA(se)) {
    ## Negative values in the Hessian matrix suggest that the algorithm did not converge
    convergence <- FALSE
    convergDiagnostic <- "nlminb optimizer reached relative convergence, saddle point?"
    warning(moduleName, "> ", convergDiagnostic, immediate. = TRUE)
  } else {
    convergDiagnostic <- out$message
  }

  ## Parameters scaling: Revert back estimated coefficients to their original scale
  out$par <- drop(out$par %*% sm)

  l <- list(formula = as.formula(fireSense_ignitionFormula),
            family = family,
            coef = setNames(out$par[1:nx], colnames(mm)),
            coef.se = setNames(se[1:nx], colnames(mm)),
            LL = -out$objective,
            AIC = 2 * length(out$par) + 2 * out$objective,
            convergence = convergence,
            convergenceDiagnostic = convergDiagnostic)

  if (hvPW) {
    l$knots <- setNames(out$par[(nx + 1L):(nx + nk)], kNames)
    l$knots.se <- setNames(se[(nx + 1L):(nx + nk)], kNames)
  }

  if (isFamilyNB) {
    ## Update the NB family template with the estimated theta
    l$family <- MASS::negative.binomial(theta = out$par[length(out$par)], link = family$link)
    l$theta <- out$par[length(out$par)]
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
