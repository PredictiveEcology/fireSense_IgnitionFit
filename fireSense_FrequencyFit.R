# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_FrequencyFit",
  description = "Fit statistical models that can be used to parametrize 
  the fire ignition component of simulation models (e.g. fireSense)",
  keywords=c("fire frequency", "optimization", "poisson", "negative binomial", "fireSense"),
  authors=c(person("Jean", "Marchal", email="jean.d.marchal@gmail.com", role=c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("0.0.0.9000"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_FrequencyFit.Rmd"),
  reqdPkgs = list("DEoptimR", "MASS", "magrittr", "numDeriv"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter(name = "formula", class = "formula", default = NULL,
      desc = 'an object of class "formula" : a symbolic description of the model to be fitted.
              Piecewised terms can be specifed using pw(variableName, knotName).'),
    defineParameter(name = "family", class = "function, character", default = "negative.binomial",
      desc = "see ?family, can be a function or a character string"),
    defineParameter(name = "start", class = "numeric, list", default = NULL, 
      desc = "optional starting values for the parameters to be estimated. Those are passed to nlminb
              and can be a numeric vector, or a list of numeric vectors. In the latter case, only the
              best solution, i.e. which minimized the objective function the most, will be kept."),
    defineParameter(name = "lb", class = "list", default = NULL, 
      desc = "optional list of numeric values describing lower bounds for the parameters to be 
              estimated. List elements are all optional but should be named (if supplied) as 'beta',
              'theta' (ignored if family is not set to negative.binomial) and 'knots'. Partial 
              matching is allowed. Values are recycled if necessary."),
    defineParameter(name = "ub", class = "numeric", default = NULL, 
      desc = "optional list of numeric values describing lower bounds for the parameters to be 
              estimated. List elements are all optional but should be named (if supplied) as 'beta',
              'theta' (ignored if family is not set to negative.binomial) and 'knots'. Partial 
              matching is allowed. Values are recycled if necessary."),
    defineParameter(name = "nlminb.control", class = "numeric", default = list(iter.max = 5e3L, eval.max=5e3L),
      desc = "optional list of control parameters to be passed to the nlminb optmizer. See ?nlminb"),
    defineParameter(name = "trace", class = "numeric", default = 0,
      desc = "non-negative integer. If > 0, tracing information on the progress of the optimization is
              produced every trace iteration. Defaults to 0 which indicates no trace information is to
              be printed."),
    defineParameter(name = "data", class = "character", default = NULL,
      desc = "optional. A character vector indicating the names of objects present in the simList 
              environment, in which to look for variables with which to predict. Objects can be
              data.frames. If omitted, or if variables are not found in data objects, variables are
              searched in the simList environment."),
    defineParameter(name = "initialRunTime", class = "numeric", default = NA, 
      desc = "optional. Simulation time at which to start this module. If omitted, start at start(simList)."),
    defineParameter(name = "intervalRunModule", class = "numeric", default = NA, 
      desc = "optional. Interval in simulation time units between two module runs.")),
  inputObjects = data.frame(objectName = "dataFireSense_FrequencyFit",
                            objectClass = "data.frame",
                            other=NA_character_,
                            stringsAsFactors=FALSE),
  outputObjects = data.frame(objectName="fireSense_FrequencyFitted",
                             objectClass="fireSense_FrequencyFit",
                             other=NA_character_,
                             stringsAsFactors=FALSE)
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_FrequencyFit = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {

    sim <- sim$fireSense_FrequencyFitInit(sim)
    
  } else if (eventType == "run") {
    
    sim <- sim$fireSense_FrequencyFitRun(sim)

  } else if (eventType == "save") {
    # ! ----- EDIT BELOW ----- ! #
    # do stuff for this event
    
    # e.g., call your custom functions/methods here
    # you can define your own methods below this `doEvent` function
    
    # schedule future event(s)
    
    # e.g.,
    # sim <- scheduleEvent(sim, time(sim) + increment, "fireSense_SizeFit", "save")
    
    # ! ----- STOP EDITING ----- ! #
    
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
fireSense_FrequencyFitInit <- function(sim) {

  sim <- scheduleEvent(sim, eventTime = if (is.na(p(sim)$initialRunTime)) start(sim) else p(sim)$initialRunTime, "fireSense_FrequencyFit", "run")
  sim

}

fireSense_FrequencyFitRun <- function(sim) {
  
  ## Toolbox: set of functions used internally by fireSense_FrequencyFitRun
    ## Handling piecewise terms in a formula
      pw <- function(v, k) pmax(v - k)
      
    ## Compute the order of magnitude
      oom <- function(x) 10^(ceiling(log10(abs(x))))
    
    ## Extract the components of the special terms, i.e. the variable and the knot value
      extractSpecial <- function(variableName, knotName) {
        cl <- match.call()

        if(missing(knotName))
          stop(paste0("fireSense_FrequencyFit> argument 'knotName' is missing (variable '", as.character(cl$v), "')"))
        else
          list(variable = as.character(cl$v), knot = as.character(cl$k))
      }
    
    ## Function to pass to the optimizer
      obj <- function(params, family, nll, scalMx, nX, mm, envData) {
        
        ## Parameters scaling
        params <- drop(params %*% scalMx)
        
        mu <- drop(mm %*% params[1L:nX])
        
        ## link implementation
        mu <- family$linkinv(mu)
        
        if(any(mu <= 0) || anyNA(mu) || any(is.infinite(mu)) || length(mu) == 0){
          return(1e20)
        } else {
          return(eval(nll, envir = as.list(environment()), enclos = envData))
        }
      }
      
    ## Function to pass to the optimizer (PW version)
      objPW <- function(params, formula, family, nll, scalMx, updateKnotExpr, nX, envData){

        ## Parameters scaling
        params <- drop(params %*% scalMx)
        
        eval(updateKnotExpr, envir = environment(), enclos = envData) ## update knot's values

        mu <- drop(model.matrix(formula, envData) %*% params[1:nX])
        
        ## link implementation
        mu <- family$linkinv(mu)

        if(any(mu <= 0) || anyNA(mu) || any(is.infinite(mu)) || length(mu) == 0){
          return(1e20)
        } else {
          return(eval(nll, envir = as.list(environment()), enclos = envData))
        }
      }
      
    ## Nlminb wrapper
      objNlminb <- function(start, objective, lower, upper, control) {
        
        nlminb.call <- quote(nlminb(start = start, objective = objective, lower = lower, upper = upper, control = control))
        nlminb.call[names(formals(objective)[-1L])] <- parse(text = formalArgs(objective)[-1L])

        o <- eval(nlminb.call)
        
        i <- 1L
        
        ## When there is no convergence and restart is possible, run nlminb() again
        while(as.integer(gsub("[\\(\\)]", "", regmatches(o$message, gregexpr("\\(.*?\\)", o$message))[[1L]])) %in% 7:14 & i < 3L){
          i <- i + 1L
          o <- eval(nlminb.call)
        }            
        o
      }
      
  envData <- new.env(parent = envir(sim))
  on.exit(rm(envData))
  list2env(as.list(envir(sim)), envir = envData)
  envData$pw <- pw
  
  if (!is.null(p(sim)$data)) ## Handling data arg
    lapply(p(sim)$data, function(x, envData) if (is.list(sim[[x]])) list2env(sim[[x]], envir = envData), envData = envData)

  if (is.empty.model(p(sim)$formula))
    stop("fireSense_FrequencyFit> The formula describes an empty model.")
  
  terms <- terms.formula(formula <- p(sim)$formula, specials = "pw")
  
  if (attr(terms, "response"))
    y <- as.character(formula[[2L]])
  else
    stop("fireSense_FrequencyFit> Incomplete formula, the LHS is missing.")
    

  family <- p(sim)$family
  
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  
  if(is.null(formals(family)$theta)) { ## Check if family is a negative.binomial
    
    family <- family()
  
  } else {
    
    family <- family(theta = NA)
    svTheta <- if (is.null(p(sim)$svTheta)) length(envData[[y]]) / sum((envData[[y]] / mean(envData[[y]]) - 1) ^ 2)
               else p(sim)$svTheta
  }
    
  
  nX <- length(labels(terms)) + attr(terms, "intercept") ## Number of variables (covariates)
  allXY <- all.vars(terms)
  
  if (is.null(attr(terms, "specials")$pw)) {
    
    allX <- allXY[allXY != y] 
    objFun <- obj
    knotNames <- kLB <- kUB <- NULL
    
  } else { ## Presence of at least one piecewise term
    
    objFun <- objPW
    
    specialsInd <- which(unlist(lapply(attr(terms,"variables"), is.call)))
    specialsCalls <- attr(terms,"variables")[specialsInd]
    
    ## Extract the names of the knots (breakpoints)
    ## Alternative way: all.vars(terms)[!all.vars(terms) %in% rownames(attr(terms,"factors"))]
    specialsTerms <- lapply(specialsCalls, function(specialsCall){
      specialsCall[[1L]] <- quote(extractSpecial)
      eval(specialsCall)
    })
    
    knotNames <- sapply(specialsTerms, "[[", "knot")
    
    if (anyDuplicated(knotNames)) stop("fireSense_FrequencyFit> Knot's names are not unique.")
    
    nKnots <- length(knotNames)
    allX <- allXY[!allXY %in% c(y, knotNames)]
    
    ## Covariates that have a breakpoint
    pwVarNames <- sapply(specialsTerms, "[[", "variable")
    
    kUB <- if (is.null(p(sim)$ub$k))
      lapply(pwVarNames, function(x) max(envData[[x]])) %>% unlist
    else
      rep_len(p(sim)$ub$k, nKnots) ## User-defined bounds (recycled if necessary)
    
    kLB <- if (is.null(p(sim)$lb$k)) {
      lapply(pwVarNames, function(x) min(envData[[x]])) %>% unlist
    } else {
      rep_len(p(sim)$lb$k, nKnots) ## User-defined bounds (recycled if necessary)
    }
    
    invisible(mapply(knotNames, z = pwVarNames, FUN = function(w, z) envData[[w]] <- mean(envData[[z]]), SIMPLIFY = FALSE))
    
    
    updateKnotExpr <- parse(text = paste0("envData[[\"", knotNames, "\"]] = params[", (nX + 1L):(nX + nKnots), "]", collapse="; "))
  }

  mm <- model.matrix(terms, envData)
  
  ## Define the scaling matrix. This is used later in the optimization process
  ##to rescale parameter values between 0 and 1, i.e. put all variables on the same scale.
    n <- nX + nKnots
    if (exists("svTheta")) n <- n + 1L
    
    scalMx <- matrix(0, n, n)
    diag(scalMx) <- 1
  

  ## Define parameter bounds automatically if they are not supplied by user
  ## First defined the bounds for DEoptim, the first optimizer    
    switch(family$link,
           log = {
             DEoptimUB <- c(
               if (is.null(p(sim)$ub$b)) {
                 ## Automatically estimate an upper boundary for each parameter       
                 (glm(formula = formula,
                      y = FALSE,
                      model = FALSE,
                      data = envData,
                      family = poisson) %>%
                    coef %>%
                    oom(.)) * 10L
               } else {
                 rep_len(p(sim)$ub$b, nX) ## User-defined bounds (recycled if necessary)
               }, kUB)
             
             DEoptimLB <- c({
               if (is.null(p(sim)$lb$b))
                 -DEoptimUB[1L:nX] ## Automatically estimate a lower boundary for each parameter 
               else 
                 rep_len(p(sim)$lb$b, nX) ## User-defined bounds (recycled)
               }, kLB)
           }, identity = {
             DEoptimUB <- c(
               if (is.null(p(sim)$ub$b)) {
                 ## Automatically estimate an upper boundary for each parameter
                 (glm(formula = formula, 
                      y = FALSE,
                      model = FALSE,
                      data = envData,
                      family = poisson(link = "identity")) %>%
                    suppressWarnings %>%
                    tryCatch(error = function(e) lm(formula,
                                                    model = FALSE,
                                                    qr = FALSE,
                                                    data = envData %>%
                                                      unlist %>% 
                                                      setNames(knotNames) %>%
                                                      as.list)) %>%
                    coef %>%
                    oom(.)) * 10L
               } else {
                 rep_len(p(sim)$ub$b) ## User-defined bounds (recycled if necessary)
               }, kUB)
             
             DEoptimLB <- c({
               if (is.null(p(sim)$lb$b))
                 rep_len(1e-30, nX) ## Ensure non-negativity
               else
                 rep_len(p(sim)$lb$b, nX) ## User-defined bounds (recycled if necessary)
               }, kLB)
           }, stop(paste("fireSense_FrequencyFit> Link function", p(sim)$link, "is not supported.")))
    
    ## If negative.binomial family needs to add bounds for theta parameter
      if (exists("svTheta")) {
        thetaEstimate <- glm.nb(formula = formula,
                                y = FALSE,
                                model = FALSE,
                                data = envData) %>%
          suppressWarnings %>%
          `[[` ("theta") %>%
          tryCatch(error = function(e)stop("fireSense_FrequencyFit> Starting value for theta cannot be estimated, please supply a starting value."))
        
        ## Add upper and lower bounds for theta
        DEoptimUB <- c(DEoptimUB, if (is.null(p(sim)$ub$t)) 2L * thetaEstimate else p(sim)$ub$t)
        DEoptimLB <- c(DEoptimLB, if (is.null(p(sim)$lb$t)) 1e-30 else p(sim)$lb$t)
      }
    
  ## Then, define lower and upper bounds for the second optimizer (nlminb)
    nlminbUB <- if (is.null(p(sim)$ub$b)) rep_len(Inf, length(DEoptimUB)) else DEoptimUB
    
    nlminbLB <- if (is.null(p(sim)$lb$b)) {
      c(switch(family$link,
               log = rep_len(-Inf, nX),        ## log-link, default: -Inf for terms and 0 for breakpoints/knots
               identity = rep_len(1e-30, nX)), ## identity link, default: enforce non-negativity
        kLB)
    } else {
      DEoptimLB ## User-defined lower bounds for parameters to be estimated
    }
  
  ## If negative.binomial family add bounds for the theta parameter
  nlminbLB <- if (exists("svTheta") && is.null(p(sim)$lb$t))
    c(nlminbLB, 1e-30) ## Enforce non-negativity
  else if (exists("svTheta"))
    c(nlminbLB, p(sim)$lb$t)

  ## Define the log-likelihood function (objective function)
  nll <- switch(family$family,
                poisson = parse(text=paste0("-sum(dpois(x=", y,", lambda = mu, log = TRUE))")),
                parse(text=paste0("-sum(dnbinom(x=", y,", mu = mu, size = params[length(params)], log = TRUE))")))
  
  trace <- if (p(sim)$trace < 0) 0 else p(sim)$trace ## No tracing if trace < 0
  
  ## If starting values are not supplied
    if (is.null(p(sim)$start)) {
      ## First optimizer, get rough estimates of the parameter values
      ## Use these estimates to compute the order of magnitude of these parameters

      JDE <- list(iter = 0L)
      i <- 0L
      while(JDE$iter == 0L && i < 30) {
        i <- i + 1L
        JDE.call <- quote(JDEoptim(fn = objFun, lower = DEoptimLB, upper = DEoptimUB, trace = if(trace > 0) TRUE else FALSE, triter = trace))
        JDE.call[names(formals(objFun)[-1])] <- parse(text = formalArgs(objFun)[-1])
        JDE <- suppressWarnings(eval(JDE.call))
      }
      
      ## Update scaling matrix
      diag(scalMx) <- oom(JDE$par)
      
      ## Update the bounds for the knots
        if (!is.null(knotNames)) {
          
            kLB <- drop(DEoptimLB %*% solve(scalMx))[(nX + 1L):(nX + nKnots)]
            nlminbLB[(nX + 1L):(nX + nKnots)] <- if (is.null(p(sim)$lb$k)) kLB else pmax(kLB, p(sim)$lb$k)
          
            kUB <- drop(DEoptimUB %*% solve(scalMx))[(nX + 1L):(nX + nKnots)]
            nlminbUB[(nX + 1L):(nX + nKnots)] <- if(is.null(p(sim)$ub$k)) kUB else pmin(kUB, p(sim)$ub$k)
  
        }
      
      ## Second optimization with nlminb()
      ## Brute-force to make models converge & select the best fit (according to AICc criterion)
        svList <- c(lapply(1:10,function(i)pmin(pmax(rnorm(length(JDE$par),0L,2L)/10 + unname(JDE$par/oom(JDE$par)), nlminbLB), nlminbUB)),
                    list(unname(JDE$par/oom(JDE$par))))
        
        out <- lapply(svList, objNlminb, objective = objFun, lower = nlminbLB, upper = nlminbUB, control = c(p(sim)$nlminb.control, list(trace = trace)))

        ## Select best minimum amongst all trials
        out <- out[[which.min(sapply(out, "[[", "objective"))]]
    
    
  } else if (is.list(p(sim)$start)) { ## If starting values are supplied as a list of vectors of starting values
    
    ## List of vectors of user-defined starting values
    out <- lapply(p(sim)$start, objNlminb, objective = objFun, lower = nlminbLB, upper = nlminbUB, control = c(p(sim)$nlminb.control, list(trace = trace)))
    
    ## Select best minimum amongst all trials
    out <- out[[which.min(sapply(out, "[[", "objective"))]]
    
  } else if (is.vector(p(sim)$start)) { ## If starting values are supplied as a vector of starting values
    
    out <- objNlminb(p(sim)$start, objFun, nlminbLB, nlminbUB, c(p(sim)$nlminb.control, list(trace = trace)))
    
  }
  
  ## Compute the standard errors around the estimates
    hess.call <- quote(numDeriv::hessian(func = objFun, x = out$par))
    hess.call[names(formals(objFun)[-1L])] <- parse(text = formalArgs(objFun)[-1L])
    hess <- eval(hess.call)
    se <- try(drop(sqrt(diag(solve(hess))) %*% scalMx), silent = TRUE)
  
  ## Negative values in the Hessian matrix suggest that the algorithm did not converge
  if(anyNA(se)) warning("fireSense_FrequencyFit> nlminb: algorithm did not converge", noBreaks. = TRUE)
  
  ## Parameters scaling: Revert back estimated coefficients to their original scale
  out$par <- drop(out$par %*% scalMx)

  sim$fireSense_FrequencyFitted <- list(formula = formula,
                                        family = family,
                                        coef = setNames(out$par[1:nX], colnames(mm)),
                                        knots = setNames(out$par[(nX + 1L):(nX + nKnots)], knotNames),
                                        se.coef = setNames(se[1:nX], colnames(mm)),
                                        se.knots = setNames(se[(nX + 1L):(nX + nKnots)], knotNames))
  
  if(exists("svTheta")){
    ## Update the NB family template with the estimated theta
    sim$fireSense_FrequencyFitted$family <- MASS::negative.binomial(theta = out$par[length(out$par)], link = family$link)
    sim$fireSense_FrequencyFitted$theta <- out$par[length(out$par)]
    sim$fireSense_FrequencyFitted$theta.se <- se[length(se)]
  }
  
  if (!is.na(p(sim)$intervalRunModule))
    sim <- scheduleEvent(sim, time(sim) + p(sim)$intervalRunModule, "fireSense_FrequencyFit", "run")
  
  sim  
}
