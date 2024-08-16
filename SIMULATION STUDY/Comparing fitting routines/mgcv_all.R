## Function to run all scampr models (IDM, PA and PO only) with simulated data

mgcv_all <- function(structured_data, unstructured_data, quad, pred, domain.data, prune.n = 4){
  
  ## wrapper functions to optimise for spatial range parameter #################
  gam_pa <- function(formula, data, weights = NULL, coord.names = c("x", "y"), k = 100, range.interval, opt.tolerance = 3,
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
    call.list$family <- binomial(link = "cloglog")
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
  
  # image(xtabs(lambda ~ x + y, domain.data), asp = 1, axes = F)
  # text(x = domain.data$x[domain.data$count > 0] / 100, y = domain.data$y[domain.data$count > 0] / 100, labels = domain.data$count[domain.data$count > 0])
  
  # combine the datasets for gfam()
  d2 <- domain.data[ , c("x", "y", "env", "count")]
  colnames(d2)[4] <- "resp"
  d1 <- structured_data[ , c("x", "y", "env", "present")]
  colnames(d1)[4] <- "resp"
  d1$source <- 1
  d2$source <- 2
  idat <- rbind(d1, d2)
  idat$source_f <- factor(idat$source)
  idat$po_id <- as.numeric(idat$source == 2)
  
  # adjust the prediction dataframe to include source and po_id
  pred$source <- 1
  pred$source_f <- 1
  pred$po_id <- 0
  
  # bfs <- simple_basis(10, data = domain.data)
  # quad$present <- 0
  # dat.scampr <- rbind(unstructured_data, quad)
  
  # PA only model ##############################################################
  
  # fit the model
  pa.time <- system.time(assign("pa", gam_pa(resp ~ env, data = d1, range.interval = c(min_pp_dist, max_pp_dist))))
  # pa_scampr <- scampr(formula = present ~ env, data = structured_data, include.sre = T, basis.functions = bfs, sre.approx = "laplace", model.type = "PA")
  
  # predict the mean abundance rate of the prediction points
  pa_pred.time <- system.time(assign("pa.pred", predict(pa, newdata = pred)))
  # collate the timings
  pa.times <- as.numeric(pa.time[3] + pa_pred.time[3])
  # get the estimated range
  rho_pa <- pa$smooth[[1]]$gp.defn[2]
  
  # PO only model ##############################################################

  # fit the model
  po.time <- system.time(assign("po", gam_po(resp ~ env, data = d1, range.interval = c(min_pp_dist, max_pp_dist))))
  # po_scampr <- scampr(formula = present ~ env, data = dat.scampr, include.sre = T, basis.functions = bfs, sre.approx = "laplace", model.type = "PO")
  
  # predict the mean abundance rate of the prediction points
  po_pred.time <- system.time(assign("po.pred", predict(po, newdata = pred)))
  # collate the timings
  po.times <- as.numeric(po.time[3] + po_pred.time[3])
  # get the estimated range
  rho_po <- po$smooth[[1]]$gp.defn[2]
  
  # IDM ########################################################################

  # fit the model - NOTE APPEARS THAT IF pho IS TOO LARGE THEN INTERCEPT BECOMES UNIDENTIFIABLE!
  idm.time <- system.time(assign("idm", gam(cbind(resp, source) ~ source_f + env + s(x, y, bs = "gp", k = 100, m = c(3, rho_pa)) + s(x, y, bs = "gp", k = 100, by = po_id, m = c(3, rho_po)),
                family=gfam(list(binomial(link = "cloglog"), poisson)),
                data=idat, method = "REML")
  ))
  # idm_scampr <-  scampr(present ~ env, data = dat.scampr, bias.formula = ~ 1, pa.data = structured_data, include.sre = T, basis.functions = bfs, sre.approx = "laplace", model.type = "IDM", latent.po.biasing = T, po.biasing.basis.functions = bfs)
  # predict the mean abundance rate of the prediction points
  idm_pred.time <- system.time(assign("idm.pred", predict(idm, newdata = pred)))
  # collate the timings
  idm.times <- as.numeric(idm.time[3] + idm_pred.time[3])
  
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
  ret_obj <- data.frame(MODEL = c("PA", "PO", "IDM"), FIT = rep("MGCV", 3),
                        KL = c(calc_KLdiv(pa.pred), calc_KLdiv(po.pred), calc_KLdiv(idm.pred)),
                        MAE = c(calc_MAE(pa.pred), calc_MAE(po.pred), calc_MAE(idm.pred)),
                        TIME = c(pa$timing, po$timing, idm.time[3] + pa$timing + po$timing),
                        TIME_FINAL_FIT = c(pa$timing_final_fit, po$timing_final_fit, idm.time[3]),
                        TIME_PRED = c(pa_pred.time[3], po_pred.time[3], idm_pred.time[3]),
                        RHO_PA = c(rho_pa, NA, rho_pa),
                        RHO_PO = c(NA, rho_po, rho_po)
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
  