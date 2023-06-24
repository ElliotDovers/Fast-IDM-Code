## Function to run all scampr models (IDM, PA and PO only) with simulated data

scampr_all <- function(structured_data, unstructured_data, quad, pred, domain.data, prune.n = 4){
  
  # add a presence identifier to the quadrature
  quad$present <- 0
  # created stacked unstructured_data and quad to be used as the data for scampr models
  dat.scampr <- rbind(unstructured_data, quad)

  # PA only model ##############################################################
  
  # fit the base model without SRE
  pa0 <- scampr(formula = present ~ env, data = structured_data, include.sre = F, sre.approx = "laplace", model.type = "PA")
  # perform the basis opt.
  pa <- basis.search.pa(pa0, domain.data = domain.data, return.model = T)
  # predict the mean abundance rate of the prediction points
  pa_pred.time <- system.time(assign("pa.pred", predict(pa, newdata = pred)))
  # collate the timings
  pa.times <- as.numeric(pa$cpu["basis.search"] + pa_pred.time[3])
  
  # PO only model ##############################################################
  
  # fit the base model without SRE
  po0 <- scampr(formula = present ~ env, data = dat.scampr, include.sre = F, sre.approx = "laplace", model.type = "PO")
  # perform the basis opt.
  po <- basis.search.po(po0, domain.data = domain.data, return.model = T) # VA approx is fine here as we are going for speed!
  # predict the mean abundance rate of the prediction points
  po_pred.time <- system.time(assign("po.pred", predict(po, newdata = pred)))
  # collate the timings
  po.times <- as.numeric(po$cpu["basis.search"] + po_pred.time[3])
  
  # IDM ########################################################################
  
  # # fit the base model without SRE
  # idm0 <- scampr(present ~ env, data = dat.scampr, bias.formula = ~ 1, IDM.presence.absence.df = structured_data, include.sre = F, sre.approx = "laplace", model.type = "IDM", latent.po.biasing = T)
  # # perform the basis opt.
  # idm <- do.call("basis.search.idm", list(object = idm0, domain.data = domain.data, return.model = T))
  # # predict the mean abundance rate of the prediction points
  # idm_pred.time <- system.time(assign("idm.pred", predict(idm, newdata = pred)))
  # # collate the timings
  # idm.times <- as.numeric(idm$cpu["basis.search"] + idm_pred.time[3])

  # IDM with fixed basis functions from the previous ###########################

  # # fit the base model without SRE
  # idm0 <- scampr(present ~ env, data = dat.scampr, bias.formula = ~ 1, IDM.presence.absence.df = structured_data, include.sre = F, sre.approx = "laplace", model.type = "IDM", latent.po.biasing = T)
  # # perform the basis opt.
  # if (is.null(pa$basis.functions)) { # in case the PA model has no SRE
  #   idm2 <- do.call("update", list(object = idm0, include.sre = T, basis.functions = attr(pa, "bfs")[[1]], po.biasing.basis.functions = po$basis.functions))
  # } else {
  #   idm2 <- do.call("update", list(object = idm0, include.sre = T, basis.functions = pa$basis.functions, po.biasing.basis.functions = po$basis.functions))
  # }
  # # predict the mean abundance rate of the prediction points
  # idm2_pred.time <- system.time(assign("idm2.pred", predict(idm2, newdata = pred)))
  # # collate the timings
  # idm2.times <- as.numeric(pa$cpu["basis.search"] + po$cpu["basis.search"] + idm0$cpu["opt"] + idm2$cpu["opt"] + idm2_pred.time[3])
  
  # fit the base model without SRE
  idm0 <- scampr(present ~ env, data = dat.scampr, bias.formula = ~ 1, IDM.presence.absence.df = structured_data, include.sre = F, sre.approx = "laplace", model.type = "IDM", latent.po.biasing = T)
  # perform the basis opt.
  if (is.null(pa$basis.functions)) { # in case the PA model has no SRE (we need to add in some basis functions so do the smallest)
    if (is.null(po$basis.functions)) {
      idm <- do.call("update", list(object = idm0, include.sre = T, basis.functions = attr(pa, "bfs")[[2]], latent.po.biasing = F))
    } else {
      idm <- do.call("update", list(object = idm0, include.sre = T, basis.functions = attr(pa, "bfs")[[2]], po.biasing.basis.functions = po$basis.functions))
    }
  } else {
    if (is.null(po$basis.functions)) {
      idm <- do.call("update", list(object = idm0, include.sre = T, basis.functions = pa$basis.functions, latent.po.biasing = F))
    } else {
      idm <- do.call("update", list(object = idm0, include.sre = T, basis.functions = pa$basis.functions, po.biasing.basis.functions = po$basis.functions))
    }
  }
  # predict the mean abundance rate of the prediction points
  idm_pred.time <- system.time(assign("idm.pred", predict(idm, newdata = pred)))
  # collate the timings
  idm.times <- as.numeric(pa$cpu["basis.search"] + po$cpu["basis.search"] + idm0$cpu["opt"] + idm$cpu["opt"] + idm_pred.time[3])
  
  ##############################################################################
  
  # calculate metrics and return result
  calc_KLdiv <- function(scampr.pred) {
    m.prd <- exp(scampr.pred)
    return(as.numeric(pred$quad.size %*% (pred$mu * log(pred$mu / m.prd))) - as.numeric(pred$quad.size %*% (pred$mu - m.prd)))
  }
  calc_MAE <- function(scampr.pred) {
    m.prd <- exp(scampr.pred)
    return(mean(abs((m.prd - mean(m.prd)) - (pred$mu - mean(pred$mu)))))
  }
  # create return object
  # ret_obj <- data.frame(MODEL = c("PA", "PO", "IDM", "IDM2"), FIT = rep("SCAMPR", 4),
  #            KL = c(calc_KLdiv(pa.pred), calc_KLdiv(po.pred), calc_KLdiv(idm.pred), calc_KLdiv(idm2.pred)),
  #            MAE = c(calc_MAE(pa.pred), calc_MAE(po.pred), calc_MAE(idm.pred), calc_MAE(idm2.pred)),
  #            TIME = c(sum(pa.times), sum(po.times), sum(idm.times), sum(idm2.times))
  # )
  ret_obj <- data.frame(MODEL = c("PA", "PO", "IDM"), FIT = rep("SCAMPR", 3),
                        KL = c(calc_KLdiv(pa.pred), calc_KLdiv(po.pred), calc_KLdiv(idm.pred)),
                        MAE = c(calc_MAE(pa.pred), calc_MAE(po.pred), calc_MAE(idm.pred)),
                        TIME = c(pa$cpu["opt"], po$cpu["opt"], idm$cpu["opt"]),
                        ALL_TIME = c(sum(pa.times), sum(po.times), sum(idm.times))
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
  