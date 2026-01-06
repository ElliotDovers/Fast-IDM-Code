gam_po <- function(formula, data, weights = NULL, coord.names = c("x", "y"), k = 100, range.interval, opt.tolerance = 3,
                   subset = NULL, na.action, offset = NULL, 
                   optimizer = c("outer", "newton"), control = list(), scale = 0, 
                   select = FALSE, knots = NULL, sp = NULL, min.sp = NULL, H = NULL, 
                   gamma = 1, fit = TRUE, paraPen = NULL, G = NULL, in.out = NULL, 
                   drop.unused.levels = TRUE, drop.intercept = NULL, nei = NULL, 
                   discrete = FALSE, ...) {
  
  mc <- match.call() # gets the arguments (must be updated for LGCP as below)
  call.list <- as.list(mc)
  
  # checks
  if ((!all(coord.names %in% colnames(data)))) {
    stop(paste0("One of 'coord.names', ", paste(coord.names, collapse = " or "), ", not found 'data'"))
  }
  
  # get the response variable out of the formula
  resp <- all.vars(formula[[2]])
  data$response <- data[ , resp]
  # evaluate the k value within the function so it can be deparsed
  call.list$k <- k
  
  # alter the call
  call.list$family <- poisson()
  call.list$data <- data
  call.list$method <- "REML"
  # update the formula for an initial fit
  call.list$formula <- as.formula(paste0("response ~ ", as.character(formula)[3], " + s(", paste(coord.names, collapse = ", "), ", bs=\"gp\", k=", deparse(call.list$k), ", m=3)"))
  # remove the function of the call
  call.list[[1]] <- NULL
  call.list$k <- NULL
  # fit an initial model to obtain warm starting parameters
  time0 <- system.time(assign("init.mod", do.call(mgcv::gam, call.list)))
  # return(list(init.mod, time0))
  
  # re-assign k
  k <- init.mod$smooth[[length(init.mod$smooth)]]$bs.dim
  
  warm.starts <- init.mod$coefficients
  warm.starts[(length(init.mod$coefficients) - (k - 1) + 1):length(init.mod$coefficients)] <- 0 # retain only the non-smooth coefficients
  # update the starting parameters
  call.list$start <- warm.starts
  # set up the object function to be optimized
  objective_fn <- function(rho) {
    call.list$formula <- as.formula(paste0("response ~ ", as.character(formula)[3], " + s(", paste(coord.names, collapse = ", "), ", bs=\"gp\", k=", deparse(k), ", m=c(3,", deparse(rho), "))"))
    tmp.m = do.call(mgcv::gam, call.list)
    return(tmp.m$gcv.ubre) # the "method" specific criterion
  }
  # calculate the optimum
  time1 <- system.time(assign("opt", optimize(objective_fn, interval = range.interval, tol = opt.tolerance)))
  # adjust the formula for the optimized spatial range parameter
  call.list$formula <- as.formula(paste0("response ~ ", as.character(formula)[3], " + s(", paste(coord.names, collapse = ", "), ", bs=\"gp\", k=", deparse(k), ", m=c(3,", deparse(opt$minimum), "))"))
  # fit the final model
  time2 <- system.time(assign("res", do.call(mgcv::gam, call.list)))
  res$timing <- time0[3] + time1[3]
  res$timing_final_fit <- time2[3]
  return(res)
}
##############################################################################