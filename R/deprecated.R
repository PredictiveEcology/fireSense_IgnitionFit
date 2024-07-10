hockeyStickApproach <- function(sim, terms, fireSense_ignitionCovariates,
                        fireSense_ignitionFormula) {
  .Deprecated(...)

  lb <- P(sim)$lb
  ub <- P(sim)$ub

  hvPW <- !is.null(attr(terms, "specials")$pw)

  ## assign default
  knotTerm <- NULL

  if (hvPW) {
    # if lb and ub don't have knots, and hvPW is true, add bounds
    # this is necessary as of 2023 given that preceding dataPrep can be run in tandem with ignitionFit.
    # As a consequence, users may not know the range of climate variables (e.g. MDC) in advance
    # for simplicity we assign default bounds to the value of the knots as the 5% and 80% percentile.
    # This situation is much more likely than one where no bounds are desired
    # bounds must be assigned before rescaling

    specialsInd <- which(unlist(lapply(attr(terms,"variables"), is.call)))
    specialsCalls <- attr(terms,"variables")[specialsInd]

    ## Extract the names of the knots (breakpoints)
    ## Alternative way: all.vars(terms)[!all.vars(terms) %in% rownames(attr(terms,"factors"))]
    specialsTerms <- lapply(specialsCalls, function(specialsCall) {
      if (specialsCall[[1L]] == "pw") {
        specialsCall[[1L]] <- quote(extractSpecial)
        eval(specialsCall)
      }
    })

    specialsTerms <- specialsTerms[!unlist(lapply(specialsTerms, is.null))]

    ## save covariate original ranges first
    ## extract variable names
    specialVars <- rownames(attr(terms, "factors"))[attr(terms, "specials")$pw]
    notSpecialVars <- colnames(attr(terms, "factor"))
    for (x in specialVars) {
      notSpecialVars <- sub(x, "", notSpecialVars, fixed = TRUE)
    }
    notSpecialVars <- unique(unlist(strsplit(notSpecialVars, ":")))

    knotTerm <- unique(sapply(specialsTerms, "[[", "variable"))
  }
  ## this assigns a default coefficient boundary of 20 - this may be scale dependent...
  ub <- checkForNullBounds(ub, 20, 0.80, knot = knotTerm, data = fireSense_ignitionCovariates)
  lb <- checkForNullBounds(lb, 0, 0.05, knot = knotTerm, data = fireSense_ignitionCovariates)

  sim$covMinMax_ignition <- fireSense_ignitionCovariates[, lapply(.SD, range), .SDcols = notSpecialVars]

  ## check for NAs
  if (any(is.na(sim$covMinMax_ignition))) {
    stop("There are NAs in fireSense_ignitionCovariates' variables used for model. Please remove NAs.")
  }

  ## rescale variable and knots.
  if (isTRUE(P(sim)$rescaleVars)) {
    if (all(is.na(P(sim)$rescalers))) {
      ## TODO: lapply through each element in rescalers and assess which elements are to be rescaled vs normalized
      message("Variables outside of [0,1] range will be rescaled to [0,1]")

      needRescale <- fireSense_ignitionCovariates[, vapply(.SD, FUN = function(x) all(inrange(na.omit(x), 0, 1)),
                                                           FUN.VALUE = logical(1)),
                                                  .SDcols = notSpecialVars]
      needRescale <- names(needRescale)[which(!needRescale)]

      message(paste("rescaling", needRescale))
      ## knots need to be added for rescaling
      rescaleDT <- fireSense_ignitionCovariates[, ..needRescale]
      knotsDT <- rbind(as.data.table(lb$knots), as.data.table(ub$knots), idcol = "knotType")
      if (sum(dim(knotsDT)) ==  0) {
        vars <- c(needRescale, "knotType")
        knotsDT[, (vars) := numeric(0)]
      } else {
        knotsDT[, knotType := ifelse(knotType == 1, "lb", "ub")]
      }

      rescaleDT <- rbind(rescaleDT, knotsDT, fill = TRUE)

      rescaleDT[, (needRescale) := lapply(.SD, FUN = function(x) {
        fireSenseUtils::rescale(x, to = c(0, 1))
      }), .SDcols = needRescale]

      fireSense_ignitionCovariates[, (needRescale) := rescaleDT[is.na(knotType), .SD, .SDcols = needRescale]]

      if (nrow(knotsDT)) {
        lb$knots <- as.list(rescaleDT[knotType == "lb", ..needRescale])
        ub$knots <- as.list(rescaleDT[knotType == "ub", ..needRescale])
      }
    } else {
      cols <- names(P(sim)$rescalers)
      fireSense_ignitionCovariates[, (cols) := mapply(FUN = function(x, vec) {x / vec},
                                                      x = .SD, vec = P(sim)$rescalers,
                                                      SIMPLIFY = FALSE),
                                   .SDcols = cols]

      # Redo parameter bounds after rescale
      lb$knots <- sapply(names(P(sim)$rescalers), FUN = function(x, knots, vec) {
        knots[[x]] / vec[x]
      }, knots = lb$knots, vec = P(sim)$rescalers, simplify = FALSE, USE.NAMES = TRUE)

      ub$knots <- sapply(names(P(sim)$rescalers), FUN = function(x, knots, vec) {
        knots[[x]] / vec[x]
      }, knots = ub$knots, vec = P(sim)$rescalers, simplify = FALSE, USE.NAMES = TRUE)
    }
  }

  ## save for later
  mod$rescales <- if (isTRUE(P(sim)$rescaleVars)) {
    if (all(is.na(P(sim)$rescalers))) { ## TODO: allow list of rescalers to be passed with mix of NA and other vals
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

  nx <- length(labels(terms)) + attr(terms, "intercept") ## Number of variables (covariates)
  allxy <- all.vars(terms)

  kLB <- kUB <- NULL

  if (hvPW) {

    objfun <- fireSenseUtils::.objFunIgnitionPW

    kNames <- sapply(specialsTerms, "[[", "knot")

    if (anyDuplicated(kNames))
      stop(moduleName, "> Knot's names are not unique.")

    nk <- length(kNames)
    allx <- allxy[!allxy %in% c(y, kNames)]

    missing <- !allxy[!allxy %in% kNames] %in% ls(fireSense_ignitionCovariates, all.names = TRUE)

    if (s <- sum(missing)) {
      stop(moduleName, "> '", allxy[!allxy %in% kNames][missing][1L], "'",
           if (s > 1) paste0(" (and ", s - 1L, " other", if (s > 2) "s", ")"),
           " not found in data objects nor in the simList environment.")
    }

    ## Covariates that have a breakpoint
    pwVarNames <- sapply(specialsTerms, "[[", "variable", simplify = FALSE)

    kUB <- if (is.null(ub$knots)) {
      lapply(pwVarNames, function(x) max(if (is(x, "AsIs")) x else fireSense_ignitionCovariates[[x]])) %>% unlist()
    } else {
      if (is.list(ub$knots) & !is.null(names(ub$knots))) {
        ## TODO: Ceres: this needs to be tested/revised if different knots are supplied for the SAME variable.
        ## TODO: Ceres: clarify doc. only one knot per variable as of now.
        lapply(pwVarNames, function(x) ub$knots[[x]]) %>% unlist()
      } else {
        rep_len(ub$knots, nk) ## User-defined bounds (recycled if necessary)
      }
    }

    kLB <- if (is.null(lb$knots)) {
      lapply(pwVarNames, function(x) min(if (is(x, "AsIs")) x else fireSense_ignitionCovariates[[x]])) %>% unlist()
    } else {
      if (is.list(lb$knots) & !is.null(names(lb$knots))) {
        ## TODO: Ceres: this needs to be tested/revised if different  knots are supplied for the same variable.
        ## TODO: Ceres: clarify doc. only one knot per variable as of now.
        lapply(pwVarNames, function(x) lb$knots[[x]]) %>% unlist()
      } else {
        rep_len(lb$knots, nk) ## User-defined bounds (recycled if necessary)
      }
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

    updateKnotExpr <- parse(text = paste0("mod_env[[\"", kNames, "\"]] = params[",
                                          (nx + 1L):(nx + nk), "]", collapse = "; "))
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

  ## TODO: Ceres: using a Poisson needs more testing so it's prevented for now
  if (grepl(x = tolower(family$family), pattern = "poisson")) {
    stop("The Poisson implementation has not been thoroughly tested and is not ready to use.
         Please use default P(sim)$family")
  }

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
      ) %>% coef() %>% oom(.)) * 10L -> ub2

      if (anyNA(ub2)) {
        stop(moduleName, "> Automated estimation of upper bounds (coefs) failed, ",
             "please set the 'coef' element of the 'ub' parameter.")
      } else {
        ub2
      }
    } else {
      ## TODO: Ceres: potentially should also accommodate different coefs
      #for different variables supplied in a list.
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
               ## TODO: Ceres: potentially should also accomodate different coefs
               # for different variables supplied in a list.
               rep_len(lb[["coef"]], nx) ## User-defined bounds (recycled if necessary)
             }

           }, identity = {
             if (is.null(lb[["coef"]])) {
               rep_len(1e-16, nx) ## Ensure non-negativity
             } else {
               ## TODO: Ceres: potentially should also accomodate different coefs
               #for different variables supplied in a list.
               rep_len(lb[["coef"]], nx) ## User-defined bounds (recycled if necessary)
             }
           }, stop(moduleName, "> Link function ", family$link, " is not supported."))
  }, kLB)

  ## If negative.binomial family needs to add bounds for theta parameter
  if (isFamilyNB) {
    DEoptimUB <- c(DEoptimUB, if (is.null(ub$theta)) 2L * theta else ub$theta)
    DEoptimLB <- c(DEoptimLB, if (is.null(lb$theta)) 1e-16 else lb$theta) ## Enfore non-negativity
  }

  ## Then, define lower and upper bounds for the second optimizer (nlminb)
  ## Upper bounds
  nlminbUB <- DEoptimUB
  if (is.null(ub[["coef"]])) {
    nlminbUB[1:nx] <- rep_len(Inf, nx)
  }

  ## Lower bounds
  nlminbLB <- if (is.null(lb[["coef"]])) {
    c(switch(family$link,
             log = rep_len(-Inf, nx),       ## log-link, default: -Inf for terms and 0 for breakpoints/knots
             identity = rep_len(1e-16, nx)) ## identity link, default: enforce non-negativity
      , kLB)

    ## when not using DEoptimLB, theta bound must be added
    if (isFamilyNB) {
      nlminbLB <- c(nlminbLB, if (is.null(lb$t)) 1e-16 else lb$t)
    }
  } else {
    DEoptimLB ## User-defined lower bounds for parameters to be estimated
  }

  ## Define the log-likelihood function (objective function)
  nll <- switch(family$family,
                poisson = parse(text = paste0("-sum(dpois(x=", y, ", lambda = mu, log = TRUE))")),
                parse(text = paste0("-sum(dnbinom(x=", y, ", mu = mu, size = params[length(params)], log = TRUE))")))

  trace <- P(sim)$trace
  message("Creating cluster")
  if (P(sim)$cores > 1) {
    if (.Platform$OS.type == "unix") {
      mkCluster <- parallel::makeForkCluster
      cl <- mkCluster(P(sim)$cores)
    } else {
      mkCluster <- parallelly::makeClusterPSOCK
      cl <- mkCluster(P(sim)$cores, rscript_libs = .libPaths())
      # mkCluster <- parallel::makePSOCKcluster ## TODO: this attaches `snow` and breaks the module
      ## see warning: https://www.rdocumentation.org/packages/secr/versions/4.3.3/topics/Parallel
    }

    # cl <- mkCluster(P(sim)$cores) #parallelly requires lib path to be set when cluster is created
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
      library("MASS")

      ## limit spawning of additional threads from workers
      data.table::setDTthreads(1)
      RhpcBLASctl::blas_set_num_threads(1)
      RhpcBLASctl::omp_set_num_threads(1)
    })
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
      DEout <- Cache(DEoptim, fn = objfun, lower = DEoptimLB, upper = DEoptimUB,
                     control = do.call("DEoptim.control", control),
                     formula = sim$fireSense_ignitionFormula,
                     mod_env = fireSense_ignitionCovariates,
                     linkinv = linkinv, nll = nll, sm = sm, nx = nx,
                     offset = offset, updateKnotExpr = updateKnotExpr,
                     userTags = c(currentModule(sim), "DEoptim"),
                     omitArgs = c("userTags"))
      DEoptimBestMem <- DEout %>%  `[[`("optim") %>% `[[`("bestmem")
    } else {
      DEout <- Cache(DEoptim, fn = objfun, lower = DEoptimLB, upper = DEoptimUB,
                     control = do.call("DEoptim.control", control),
                     mod_env = fireSense_ignitionCovariates,
                     linkinv = linkinv, nll = nll, sm = sm, nx = nx, mm = mm,
                     offset = offset,
                     userTags = c(currentModule(sim), "DEoptim"),
                     omitArgs = c("userTags"))
      DEoptimBestMem <- DEout %>% `[[`("optim") %>% `[[`("bestmem")
    }

    if (FALSE) { # THIS IS ELIOT BLOCKING THIS -- I DON"T THINK IT IS NECESSARY
      ## Update scaling matrix
      diag(sm) <- oom(DEoptimBestMem)

      #TODO: I believe this is a mistake, and the length(nlminUB) subset should be replaced with
      ## Update of the lower and upper bounds of the coefficients based on the scaling matrix
      nlminbLB <- nlminbLB / diag(sm)
      nlminbUB <- nlminbUB / diag(sm)

      ## Update of the lower and upper bounds for the knots based on the scaling matrix
      ## TODO: fix the partial matching

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
      pop <- if (exists("cl", inherits = FALSE)) length(cl) else 1
      if (pop < 20) pop <- 20
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
        parallel::clusterExport(cl, c("outputPath", "basePattern"), envir = environment())

        parallel::clusterEvalQ(
          cl,
          sink(file.path(outputPath, paste0(basePattern, ".", Sys.getpid())))
        )
        pids <- unlist(parallel::clusterEvalQ(cl, Sys.getpid()))
        message("These are the pids of the spawned cores: ")
        dput(pids)
      }
      message("Starting nlminb ... ")

      if (hvPW) {
        out <- Cache(parallel::clusterApplyLB, cl = cl, x = start, fun = objNlminb, objective = objfun,
                     lower = nlminbLB, upper = nlminbUB, hvPW = hvPW,
                     linkinv = linkinv, nll = nll, sm = sm, nx = nx, mm = mm, #TODO mm may not be required with PW...
                     mod_env = fireSense_ignitionCovariates, offset = offset,
                     formula = sim$fireSense_ignitionFormula,
                     updateKnotExpr = updateKnotExpr,# cacheId = "e016b5d728ed2b6a",
                     control = c(P(sim)$nlminb.control, list(trace = trace)),
                     userTags = c(currentModule(sim), "objNlminb"),
                     omitArgs = c("x", "userTags"),
                     .functionName = "objNlminb") # don't need to know the random sample... the mm is enough
      } else {
        out <- Cache(parallel::clusterApplyLB, cl = cl, x = start, fun = objNlminb, objective = objfun,
                     lower = nlminbLB, upper = nlminbUB, hvPW = hvPW,
                     linkinv = linkinv, nll = nll, sm = sm, nx = nx, mm = mm, ## TODO mm may not be required with PW...
                     mod_env = fireSense_ignitionCovariates, offset = offset,
                     control = c(P(sim)$nlminb.control, list(trace = trace)),
                     userTags = c(currentModule(sim), "objNlminb_without_hvPW"),
                     omitArgs = c("x", "userTags"),
                     .functionName = "objNlminb") # don't need to know the random sample... the mm is enough
      }

      if (FALSE) { # THIS SECTION ALLOWS MANUAL READING OF LOG FILES
        #  MUST MANUALLY IDENTIFY THE PIDS
        if (exists("pids")) {
          aa <- lapply(paste0("~/GitHub/WBI_forecasts/outputs/AB/fireSense_IgnitionFit_spades189_20210308_trace.", pids), readLines)
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

      if (trace) parallel::clusterEvalQ(cl, sink())
    } else {
      ## see commit ce3a8f41823583f848793f2743281f643d29ef05 if error occurs - it may be benign
      if (hvPW) {
        out <- Cache(lapply, start, objNlminb, objective = objfun,
                     lower = nlminbLB, upper = nlminbUB, hvPW = hvPW,
                     linkinv = linkinv, nll = nll, sm = sm, nx = nx,
                     mm = mm, ## TODO mm may not be required with PW
                     mod_env = fireSense_ignitionCovariates, offset = offset,
                     formula = sim$fireSense_ignitionFormula,
                     updateKnotExpr = updateKnotExpr,
                     userTags = c(currentModule(sim), "objNlminb"),
                     omitArgs = c("X", "userTags"),
                     control = c(P(sim)$nlminb.control, list(trace = min(6, trace * 3))))
      } else {
        out <- Cache(lapply, start, objNlminb, objective = objfun,
                     lower = nlminbLB, upper = nlminbUB, hvPW = hvPW,
                     linkinv = linkinv, nll = nll, sm = sm, nx = nx,
                     mm = mm, ## TODO mm may not be required with PW
                     mod_env = fireSense_ignitionCovariates, offset = offset,
                     userTags = c(currentModule(sim), "objNlminb"),
                     omitArgs = c("X", "userTags"),
                     control = c(P(sim)$nlminb.control, list(trace = min(6, trace * 3))))
      }

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
                  formula = sim$fireSense_ignitionFormula,
                  offset = offset,
                  userTags = c(currentModule(sim), "hessian"),
                  omitArgs = c("userTags"))
  } else {
    hess <- Cache(numDeriv::hessian, func = objfun, x = outBest$par,
                  mod_env = fireSense_ignitionCovariates,
                  linkinv = linkinv, nll = nll,
                  sm = sm, nx = nx, mm = mm,
                  offset = offset,
                  userTags = c(currentModule(sim), "hessian"),
                  omitArgs = c("userTags"))
  }

  solvedHess <- tryCatch(solve(hess), error = function(e) NA)
  se <- suppressWarnings(tryCatch(drop(sqrt(diag(solvedHess)) %*% sm), error = function(e) NA))
  if (all(is.na(se))) {
    seSimple <- sqrt(1/diag(hess))
    if (all(!is.infinite(seSimple))) {
      warning("The hessian could not be inverted; but the 'crude estimate' of Bolker",
              " sqrt(1/diag(hess)) can and will be used.")
      se <- seSimple
    }
  }
  message("... Done")

  if (!exists("outBest", inherits = FALSE)) {
    best <- drop(as.matrix(dd[1, -(1:2)]))
  } else {
    best <- outBest$par
    colnms <- c(attr(terms, "term.labels"))
    if (hvPW) {
      colnms <- c(colnms, unlist(lapply(updateKnotExpr, function(x) x[[2]][[3]])))
    }
    if (isFamilyNB) {
      colnms <- c(colnms, "NB_theta")
    }

    names(best) <- colnms
  }

  if (anyPlotting(P(sim)$.plots)) {
    message("Plotting has not been tested thoroughly")

    ## Next line added by Ceres to avoid hardcoding,
    ##  but not necessary with suggested solution
    colName <- if (exists("specialsTerms", inherits = FALSE)) {
      unique(rbindlist(specialsTerms)$variable)
    } else {
      tt <- terms(as.formula(sim$fireSense_ignitionFormula[-(2)], env = .GlobalEnv))
      facts <- attr(tt, "factors")
      rownames(facts)[sapply(rownames(facts), function(v) length(grep(v, attr(tt, "term.labels"))) > 1)]
    }

    ## TODO: Ceres: this is not working if using formula/data different from Ian's/Tati's
    ## suggested solution, pass the original data.frame to get the variables (and potentially the max/min) and
    ## generate new data from it.

    xCeiling <- max(sim$fireSense_ignitionCovariates[[colName]]) * 1.5 #subsets full 0-1 range
    ndLong <- pwPlotData(bestParams = best,
                         ses = se, solvedHess = solvedHess,
                         formula = sim$fireSense_ignitionFormula,
                         xColName = colName, nx = nx, offset = offset,
                         linkinv = linkinv,
                         rescaler = mod$rescales,
                         rescaleVar = P(sim)$rescaleVars,
                         covMinMax_ignition = sim$covMinMax_ignition,
                         xCeiling = xCeiling)

    ## round to avoid silly decimal errors
    resInKm2 <- round(prod(res(sim$ignitionFitRTM)) / 1e6) ## 1e6 m^2 == 1 km^2
    labelToUse <- paste("Ignition rate per", resInKm2, "km^2")
    filenameToUse <- paste0("IgnitionRatePer", resInKm2, "km2_", P(sim)$.studyAreaName)

    Plots(data = ndLong, fn = pwPlot, xColName = colName,
          ggylab = labelToUse, # types = "screen",
          origXmax = max(sim$fireSense_ignitionCovariates[[colName]]), ## if supplied, adds bar to plot
          ggTitle =  paste("fireSense_IgnitionFit:", P(sim)$.studyAreaName,
                           "(", basename(outputPath(sim)), ")"),
          filename = filenameToUse)
    ## TODO: unresolved bug in Plot triggered by spaces

    ## FITTED VS OBSERVED VALUES
    fittedVals <- predictIgnition(as.formula(fireSense_ignitionFormula, env = .GlobalEnv),
                                  fireSense_ignitionCovariates,
                                  setNames(outBest$par[1:nx], colnames(mm)),
                                  1,
                                  1,
                                  family$linkinv)

    plotData <- data.table(fireSense_ignitionCovariates)
    plotData[,  rows := 1:nrow(plotData)]
    cols <- unique(c(paste(y), xvar, "rows"))
    plotData <- plotData[, ..cols]
    plotData <- cbind(plotData, fittedVals = fittedVals)

    predDT <- rbindlist(lapply(1:100,  FUN = function(x, DT, theta) {
      rnbinomPred <- rnbinom(nrow(DT), mu = DT$fittedVals,
                             size = theta)
      n <- rep(x, nrow(DT))
      data.table(rnbinomPred = rnbinomPred, n = n, rows = DT$rows)
    }, DT = plotData, theta = outBest$par[length(outBest$par)]))

    plotData <- plotData[predDT, on = "rows"]

    plotData <- plotData[, list(obsFires = sum(eval(y), na.rm = TRUE),
                                predFires = sum(rnbinomPred, na.rm = TRUE)),
                         by = c(xvar, "n")]
    plotData[, predFires := as.integer(predFires)]
    plotData <- melt(plotData, id.var = c(xvar, "n"))

    Plots(data = plotData, fn = fittedVsObservedPlot,
          xColName = xvar,
          ggylab = "num. fires",
          ggTitle = paste(P(sim)$.studyAreaName, "fireSense_IgnitionFit: obs. vs. fitted values"),
          filename = paste0("ignition_NumFiresFitted", P(sim)$.studyAreaName))
  }

  convergence <- TRUE

  if (outBest$convergence || anyNA(se)) {
    tooClose <- 0.00001

    # This includes the pw terms on their own -- but if a pw term is not estimable, we need
    #   to remove the interaction term too
    closeToBounds <- abs(drop((best - DEoptimLB) / (DEoptimUB - DEoptimLB))) < tooClose |
      abs(best) < tooClose
    ctb <- data.table(term = names(closeToBounds), best = best,
                      upperBoundary = DEoptimUB,
                      lowerBoundary = DEoptimLB,
                      closeToBounds = closeToBounds)[closeToBounds]

    closeToBoundsOnlyCovariates <- closeToBounds[1:nx]
    # This only has terms with covariates (including pw in interaction), not pw terms on their own
    possTerms <- attr(terms, "term.labels")[1:nx]
    if (nrow(ctb)) {
      toRemove <- do.call(rbind, lapply(ctb$term, function(trm) grepl(trm, possTerms)))
      toRemove <- apply(toRemove, 2, any)
    } else {
      toRemove <- sapply(closeToBoundsOnlyCovariates, function(x) FALSE)
    }
    toRemove <- toRemove | closeToBoundsOnlyCovariates

    possTerms <- possTerms[!toRemove]
    possForm <- paste0(terms[[2]], " ~ ", paste(possTerms, collapse = " + "), " -1")

    tryRefit <- TRUE
    ## TODO: this may not be the best way for checking if the formulas are identical.
    if (identical(as.character(parse(text = sim$fireSense_ignitionFormula)),
                  as.character(parse(text = possForm)))) {
      tryRefit <- FALSE
      message("Can't find a simpler model to refit automatically. Will use the current formula:")
      message(possForm)
      message("but note that there are convergence or invertability issues preventing calculation of std errors.")
    }

    message("--------------------------------------------------")
    message("It is possible that parameters are too close to their lower boundary values (or zero). ",
            "The following are within ", tooClose*100, "% of their lower boundary and removing them ",
            "from sim$fireSense_ignitionFormula may help with convergence or invertability... e.g.")
    message("P(sim)$fireSense_ignitionFit$fireSense_ignitionFormula <- \"", possForm, "\"")
    messageDF(ctb)
    message("It may also help to use Ben Bolker's approximation: sqrt(1/diag(hess)) mentioned here:")
    message("<https://cran.r-project.org/web/packages/bbmle/vignettes/mle2.pdf>")
    message("If there are Inf values, that indicates variables to remove as they have",
            "infinite variance at the solution")

    if (!isFALSE(P(sim)$autoRefit) & tryRefit)  {
      outRL <- if (isTRUE(P(sim)$autoRefit)) {
        message("Automatically refitting with simpler model because P(sim)$autoRefit is TRUE")
        "y"
      } else if (isFALSE(P(sim)$autoRefit)) {
        "n"
      } else {
        readline("Would you like to restart this IgnitionFit event with that new formula (Y or N or interactive)? ")
      }
      if (identical(tolower(outRL), "y")) {
        sim$fireSense_ignitionFormula <- possForm
        sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "run", eventPriority = 2)
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

  ## TODO: added raster attributes may not be ideal to track number of non-NAs
  ## rationale for lambdaRescaleFactor:
  ## original fire prob is sum(n_fires)/nrow(preSampleData),
  ## the fitted one, imposed by sampling, becomes sum(n_fires)/nrow(postSampleData)
  ## so to adjust predicted values, one needs to predVals * nrow(postSampleData)/nrow(preSampleData)
  origNoPix <- attributes(sim$ignitionFitRTM)$nonNAs   ## nrow(preSampleData) in eg above
  finalNoPix <- nrow(fireSense_ignitionCovariates)     ## nrow(postSampleData) in eg above
  lambdaRescaleFactor <- finalNoPix/origNoPix

  l <- list(formula = as.formula(fireSense_ignitionFormula, env = .GlobalEnv),
            family = family,
            data = fireSense_ignitionCovariates,
            coef = setNames(outBest$par[1:nx], colnames(mm)),
            coef.se = setNames(se[1:nx], colnames(mm)),
            LL = -outBest$objective,
            AIC = 2 * length(outBest$par) + 2 * outBest$objective,
            convergence = convergence,
            convergenceDiagnostic = convergDiagnostic,
            rescales = mod$rescales,
            fittingRes = res(sim$ignitionFitRTM)[1],
            lambdaRescaleFactor = lambdaRescaleFactor)

  if (hvPW) {
    l$knots <- setNames(outBest$par[(nx + 1L):(nx + nk)], kNames)
    l$knots.se <- setNames(se[(nx + 1L):(nx + nk)], kNames)
  }

  if (isFamilyNB) {
    ## Update the NB family template with the estimated theta
    l$family <- MASS::negative.binomial(theta = outBest$par[length(outBest$par)], link = family$link)
    l$theta <- outBest$par[length(outBest$par)]
    l$theta.se <- se[length(se)]

    return(l)
  }
}
