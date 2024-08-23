## Function to run all scampr models (IDM, PA and PO only) with simulated data

mgcv_test <- function(structured_data, unstructured_data, quad, pred, domain.data, prune.n = 4){
  
  ## wrapper functions to optimise for spatial range parameter #################
  ## USING POISSON GRID COUNTS #################################################
  gam_lgcp <- function(formula, data, weights = NULL, coord.names = c("x", "y"), k = 100, range.interval, opt.tolerance = 3,
                       subset = NULL, na.action, offset = NULL, 
                       optimizer = c("outer", "newton"), control = list(), scale = 0, 
                       select = FALSE, knots = NULL, sp = NULL, min.sp = NULL, H = NULL, 
                       gamma = 1, fit = TRUE, paraPen = NULL, G = NULL, in.out = NULL, 
                       drop.unused.levels = TRUE, drop.intercept = NULL, nei = NULL, 
                       discrete = FALSE, ...) {
    
    mc <- match.call() # gets the arguments (must be updated for LGCP as below)
    call.list <- as.list(mc)
    
    # check the form of the weights
    object.supplied <- tryCatch(!is.null(weights), error = function(e) FALSE)
    if (!object.supplied) {
      weight.name <- deparse(substitute(weights))
      if (weight.name == "NULL") {
        stop("weights must be supplied for fitting a LGCP.\n\nThese must be quadrature weights in rows where the formula response = 0,\nThese must be some small number (e.g. 1e-6) in rows where the formula response = 1.")
      } else {
        wt.vec <- as.vector(data[,weight.name])
      }
    } else {
      wt.vec <- weights
    }
    # checks
    if ((!all(coord.names %in% colnames(data)))) {
      stop(paste0("One of 'coord.names', ", paste(coord.names, collapse = " or "), ", not found 'data'"))
    }
    
    # get the response variable out of the formula
    resp <- all.vars(formula[[2]])
    data$new.response <- data[ , resp] / wt.vec
    
    # alter the call according to requirements for a LGCP
    call.list$family <- poisson()
    call.list$data <- data
    call.list$weights <- wt.vec
    call.list$method <- "REML"
    # update the formula for an initial fit
    call.list$formula <- as.formula(paste0("new.response ~ ", as.character(formula)[3], " + s(", paste(coord.names, collapse = ", "), ", bs=\"gp\", k=", deparse(k), ", m=3)"))
    # remove the function of the call
    call.list[[1]] <- NULL
    # fit an initial model to obtain warm starting parameters
    init.mod <- do.call(mgcv::gam, call.list)
    warm.starts <- init.mod$coefficients
    warm.starts[(length(init.mod$coefficients) - k + 1):length(init.mod$coefficients)] <- 0
    # update the starting parameters
    call.list$start <- warm.starts
    # set up the object function to be optimized
    objective_fn <- function(rho) {
      call.list$formula <- as.formula(paste0("new.response ~ ", as.character(formula)[3], " + s(", paste(coord.names, collapse = ", "), ", bs=\"gp\", k=", deparse(k), ", m=c(3,", deparse(rho), "))"))
      tmp.m = do.call(mgcv::gam, call.list)
      return(tmp.m$gcv.ubre) # the "method" specific criterion
    }
    # calculate the optimum
    opt <- optimize(objective_fn, interval = range.interval, tol = opt.tolerance)
    # adjust the formula for the optimized spatial range parameter
    call.list$formula <- as.formula(paste0("new.response ~ ", as.character(formula)[3], " + s(", paste(coord.names, collapse = ", "), ", bs=\"gp\", k=", deparse(k), ", m=c(3,", deparse(opt$minimum), "))"))
    # fit the final model
    res <- do.call(mgcv::gam, call.list)
    return(res)
  }
  ##############################################################################
  
  ## wrapper functions to optimise for spatial range parameter #################
  ## USING BERMAN-TURNER DEVICE ################################################
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
    
    # alter the call
    call.list$family <- poisson()
    call.list$data <- data
    call.list$method <- "REML"
    # update the formula for an initial fit
    call.list$formula <- as.formula(paste0("response ~ ", as.character(formula)[3], " + s(", paste(coord.names, collapse = ", "), ", bs=\"gp\", k=", deparse(k), ", m=3)"))
    # remove the function of the call
    call.list[[1]] <- NULL
    # fit an initial model to obtain warm starting parameters
    time0 <- system.time(assign("init.mod", do.call(mgcv::gam, call.list)))
    # return(list(init.mod, time0))
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
  
  ## Set up the data
  
  library(fields)
  
  # calculate the distances from points to domain quadrature centers
  dist2quad <- rdist(unstructured_data[ , c("x", "y")], domain.data[ , c("x", "y")])
  dist_between_pts <- rdist(unstructured_data[ , c("x", "y")])
  max_pp_dist <- max(dist_between_pts)
  diag(dist_between_pts) <- max_pp_dist
  min_pp_dist <- min(dist_between_pts)
  # find indices of the closest quadrat
  quad.id <- table(apply(dist2quad, 1, which.min))
  
  # set up the new Poisson count response
  domain.data$count <- 0
  domain.data$count[as.numeric(names(quad.id))] <- as.vector(quad.id)
  
  # add a presence identifier to the quadrature
  quad$present <- 0
  # change to small weights on the presence points for Berman-Turner device
  unstructured_data$quad.size <- 1e-6
  # created stacked unstructured_data and quad to be used as the data for scampr models
  dat.scampr <- rbind(unstructured_data, quad)
  
  # PO only model ##############################################################

  po.time <- system.time(assign("po", gam_lgcp(formula = present ~ env, data = dat.scampr, weights = dat.scampr$quad.size, range.interval = c(min_pp_dist, max_pp_dist))))
  
  # fit the model
  po_grid.time <- system.time(assign("po_grid", gam_po(count ~ env, data = domain.data, range.interval = c(min_pp_dist, max_pp_dist))))
  # po_scampr <- scampr(formula = present ~ env, data = dat.scampr, include.sre = T, basis.functions = bfs, sre.approx = "laplace", model.type = "PO")
  
  # predict the mean abundance rate of the prediction points
  po_pred.time <- system.time(assign("po.pred", predict(po, newdata = pred)))
  po_grid_pred.time <- system.time(assign("po_grid.pred", predict(po_grid, newdata = pred)))
  # collate the timings
  po.times <- as.numeric(po.time[3] + po_pred.time[3])
  po_grid.times <- as.numeric(po_grid.time[3] + po_grid_pred.time[3])
  # get the estimated range
  rho_po <- po$smooth[[1]]$gp.defn[2]
  rho_po_grid <- po_grid$smooth[[1]]$gp.defn[2]

  ##############################################################################
  
  # calculate metrics and return result
  calc_KLdiv <- function(mgcv.pred) {
    m.prd <- exp(mgcv.pred)
    return(as.numeric(pred$quad.size %*% (pred$mu * log(pred$mu / m.prd))) - as.numeric(pred$quad.size %*% (pred$mu - m.prd)))
  }
  calc_MAE <- function(mgcv.pred) {
    m.prd <- exp(mgcv.pred)
    return(mean(abs((m.prd - mean(m.prd)) - (pred$mu - mean(pred$mu)))))
  }
  # create return object
  # ret_obj <- data.frame(MODEL = c("PA", "PO", "IDM", "IDM2"), FIT = rep("SCAMPR", 4),
  #            KL = c(calc_KLdiv(pa.pred), calc_KLdiv(po.pred), calc_KLdiv(idm.pred), calc_KLdiv(idm2.pred)),
  #            MAE = c(calc_MAE(pa.pred), calc_MAE(po.pred), calc_MAE(idm.pred), calc_MAE(idm2.pred)),
  #            TIME = c(sum(pa.times), sum(po.times), sum(idm.times), sum(idm2.times))
  # )
  ret_obj <- data.frame(MODEL = c("BT", "GRID"), FIT = rep("MGCV", 2),
                        KL = c(calc_KLdiv(po.pred), calc_KLdiv(po_grid.pred)),
                        MAE = c(calc_MAE(po.pred), calc_MAE(po_grid.pred)),
                        TIME = c(po$timing, po_grid$timing),
                        TIME_FINAL_FIT = c(NA, NA),
                        TIME_PRED = c(po_pred.time[3], po_grid_pred.time[3]),
                        RHO_PA = c(NA, NA),
                        RHO_PO = c(rho_po, rho_po_grid)
  )
  # # alter the PA search results to combine
  # tmp.pa <- attr(pa, "search.res")
  # tmp.pa$k_bias <- NA
  # tmp.pa$radius_bias <- NA
  # tmp.pa$MODEL = "PA"
  # tmp.pa$FIT = "SCAMPR"
  # # alter the PO search results to combine
  # tmp.po <- attr(po, "search.res")
  # tmp.po$k_bias <- NA
  # tmp.po$radius_bias <- NA
  # tmp.po$MODEL = "PO"
  # tmp.po$FIT = "SCAMPR"
  # # alter the IDM search results to combine
  # # tmp.idm <- attr(idm, "search.res")
  # # tmp.idm$MODEL = "IDM"
  # # tmp.idm$FIT = "SCAMPR"
  # tmp.idm <- data.frame(nodes = c(0, NA), k = c(0, NA), radius = c(NA, NA),
  #                       ll = c(logLik(idm0), logLik(idm)), BIC = BIC(idm0, idm)$BIC,
  #                       cpu = c(idm0$cpu["opt"], idm0$cpu["opt"]), convergence = c(idm0$convergence, idm$convergence),
  #                       k_bias = c(0, NA), radius_bias = c(0, NA), MODEL = "IDM", FIT = "SCAMPR"
  # )
  # # combine
  # attr(ret_obj, "all") <- rbind(
  #   rbind(tmp.pa, tmp.po), tmp.idm
  # )
  
  return(ret_obj)
}
  