## Function to run mgcv IDM

idm_mgcv <- function(form, dat_pa, pres, quad, domain, fld_dim = 100, bias_fld_dim = 100, po.as.grid = FALSE, estimate.range = FALSE, with.predictions = FALSE){
  
  library(mgcv)
  source("gam_pa.R")
  source("gam_po.R")
  
  ## Set up the data
  
  library(fields)
  
  # calculate the intra-point distances
  dist_between_pts <- rdist(pres[ , c("x", "y")])
  max_pp_dist <- max(dist_between_pts)
  diag(dist_between_pts) <- max_pp_dist
  min_pp_dist <- min(dist_between_pts)
  
  if (po.as.grid) {
    
    # calculate the distances from points to domain quadrature centers
    dist2quad <- rdist(pres[ , c("x", "y")], domain[ , c("x", "y")])
    # find indices of the closest quadrat
    quad.id <- table(apply(dist2quad, 1, which.min))
    
    # set up the new Poisson count response
    domain$count <- 0
    domain$count[as.numeric(names(quad.id))] <- as.vector(quad.id)
    
    # image(xtabs(lambda ~ x + y, domain.data), asp = 1, axes = F)
    # text(x = domain.data$x[domain.data$count > 0] / 100, y = domain.data$y[domain.data$count > 0] / 100, labels = domain.data$count[domain.data$count > 0])
    
    # combine the datasets for gfam()
    d2 <- domain[ , c(all.vars(form[[3]]), "x", "y", "count", "quad.size")]
    colnames(d2)[colnames(d2) == "count"] <- all.vars(form[[2]])
    
  } else {
    
    # setup the required data structure
    pres$quad.size <- 1e-6
    dat_po <- rbind(pres, quad)
    
    # combine the datasets for gfam()
    d2 <- dat_po[ , c(all.vars(form[[3]]), "x", "y", "occ", "quad.size")]
    colnames(d2)[colnames(d2) == "wt"] <- "quad.size"
    
  }
  
  d1 <- dat_pa[ , c(all.vars(form[[2]]), all.vars(form[[3]]), "x", "y", "quad.size")]
  
  d1$source <- 1
  d2$source <- 2
  idat <- rbind(d1, d2)
  idat$po_id <- as.numeric(idat$source == 2)
  
  # adjust the prediction dataframe to include source and po_id
  domain$source <- 1
  domain$po_id <- 0
  
  if (estimate.range) {
    
    ### Estimate spatial range of error term from the PA only model ############
    
    # fit the model
    pa.time <- system.time(assign("pa", gam_pa(form, data = d1, k = fld_dim, range.interval = c(min_pp_dist, max_pp_dist))))
  
    # get the estimated range
    rho_pa <- pa$smooth[[1]]$gp.defn[2]
    
    # Estimate spatial range of error term from the PO only model ##############
    
    # fit the model
    po.time <- system.time(assign("po", gam_po(form, data = d2, k = bias_fld_dim, range.interval = c(min_pp_dist, max_pp_dist))))
    
    # get the estimated range
    rho_po <- po$smooth[[1]]$gp.defn[2]
    
    # IDM ########################################################################
    
    # update the formula
    form[[2]] <- bquote(cbind(occ, source))
    # new.form <- as.formula(paste0(paste(deparse(form), collapse = ""), paste0(" + factor(source) + s(x, y, bs = 'gp', k = ", fld_dim, ", m = c(3, rho_pa)) + s(x, y, bs = 'gp', k = ", bias_fld_dim, ", by = po_id, m = c(3, rho_po))")))
    new.form <- as.formula(paste0(paste(deparse(form), collapse = ""), paste0(" + s(x, y, bs = 'gp', k = ", fld_dim, ", m = c(3, rho_pa)) + s(x, y, bs = 'gp', k = ", bias_fld_dim, ", by = po_id, m = c(3, rho_po))")))
    
  } else {
    
    # update the formula
    form[[2]] <- bquote(cbind(occ, source))
    # new.form <- as.formula(paste0(paste(deparse(form), collapse = ""), paste0(" + factor(source) + s(x, y, bs = 'gp', k = ", fld_dim, ", m = c(3, rho_pa)) + s(x, y, bs = 'gp', k = ", bias_fld_dim, ", by = po_id, m = c(3, rho_po))")))
    new.form <- as.formula(paste0(paste(deparse(form), collapse = ""), paste0(" + s(x, y, bs = 'gp', k = ", fld_dim, ", m = c(3)) + s(x, y, bs = 'gp', k = ", bias_fld_dim, ", by = po_id, m = c(3))")))
    
    
  }
  
  # fit the model
  idm.time <- system.time(assign("idm", gam(new.form, family=gfam(list(binomial(link = "cloglog"), poisson)), data=idat, offset = log(idat$quad.size), method = "REML")
  ))

  if (with.predictions) {
    # predict the mean abundance rate of the prediction points
    idm$pred <- predict(idm, newdata = domain)
  }
  
  return(idm)
}
  