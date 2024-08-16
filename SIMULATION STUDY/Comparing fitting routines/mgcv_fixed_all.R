## Function to run all scampr models (IDM, PA and PO only) with simulated data

mgcv_fixed_all <- function(structured_data, unstructured_data, quad, pred, domain.data, prune.n = 4){
  
  library(fields)
  
  # calculate the distances from points to domain quadrature centers
  dist2quad <- rdist(unstructured_data[ , c("x", "y")], domain.data[ , c("x", "y")])
  max_pp_dist <- max(rdist(unstructured_data[ , c("x", "y")]))
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
  
  # PA only model ##############################################################
  
  # fit the model
  pa.time <- system.time(assign("pa", gam(resp ~ env + s(x, y, bs = "gp", k = 100, m = c(3, max_pp_dist)),
                                            family=binomial(link = "cloglog"),
                                            data=d1, method = "REML")
  ))
  # pa_scampr <- scampr(formula = present ~ env, data = structured_data, include.sre = T, basis.functions = bfs, sre.approx = "laplace", model.type = "PA")
  
  # predict the mean abundance rate of the prediction points
  pa_pred.time <- system.time(assign("pa.pred", predict(pa, newdata = pred)))
  # collate the timings
  pa.times <- as.numeric(pa.time[3] + pa_pred.time[3])
  
  # PO only model ##############################################################

  # fit the model
  po.time <- system.time(assign("po", gam(resp ~ env + s(x, y, bs = "gp", k = 100, m = c(3, max_pp_dist)),
                                          family=poisson,
                                          data=d2, method = "REML")
  ))
  # po_scampr <- scampr(formula = present ~ env, data = dat.scampr, include.sre = T, basis.functions = bfs, sre.approx = "laplace", model.type = "PO")
  # predict the mean abundance rate of the prediction points
  po_pred.time <- system.time(assign("po.pred", predict(po, newdata = pred)))
  # collate the timings
  po.times <- as.numeric(po.time[3] + po_pred.time[3])
  
  # IDM ########################################################################

  # fit the model - NOTE APPEARS THAT IF pho IS TOO LARGE THEN INTERCEPT BECOMES UNIDENTIFIABLE!
  idm.time <- system.time(assign("idm", gam(cbind(resp, source) ~ source_f + env + s(x, y, bs = "gp", k = 100, m = c(3, 50)) + s(x, y, bs = "gp", k = 100, by = po_id, m = c(3, 50)),
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
                        TIME = c(pa.time[3], pa.time[3], idm.time[3]),
                        TIME_FINAL_FIT = c(pa.time[3], po.time[3], idm.time[3]),
                        TIME_PRED = c(pa_pred.time[3], po_pred.time[3], idm_pred.time[3]),
                        RHO_PA = c(max_pp_dist, NA, max_pp_dist),
                        RHO_PO = c(NA, max_pp_dist, max_pp_dist)
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
  