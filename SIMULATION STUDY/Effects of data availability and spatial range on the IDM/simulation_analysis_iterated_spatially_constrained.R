home.wd <- getwd()

## Install required packages ###################################################
if(!require(scampr, quietly = T)){
  # scampr package can be installed from source provided in the code zip
  setwd("..")
  setwd("..")
  install.packages(paste0(getwd(), "/scampr_0.0.0.9000.tar.gz"), repos = NULL, type="source")
  library(scampr)
  setwd(home.wd)
}
if(!require(RandomFieldsUtils, quietly = T)){
  # RandomFieldsUtils package can be installed from source provided in the code zip
  setwd("..")
  setwd("..")
  install.packages(paste0(getwd(), "/RandomFieldsUtils_1.2.5.tar.gz"), repos = NULL, type="source")
  library(RandomFieldsUtils)
  setwd(home.wd)
}
if(!require(RandomFields, quietly = T)){
  # RandomFields package can be installed from source provided in the code zip
  setwd("..")
  setwd("..")
  install.packages(paste0(getwd(), "/RandomFields_3.3.14.tar.gz"), repos = NULL, type="source")
  library(RandomFields)
  setwd(home.wd)
}
if(!require(sp, quietly = T)){
  install.packages("sp")
  library(sp)
}
if(!require(fields, quietly = T)){
  install.packages("fields")
  library(fields)
}
if(!require(spatstat, quietly = T)){
  install.packages("spatstat")
  library(spatstat)
}
################################################################################

# Get the job array
tab <- read.csv("job_array.csv")

# TOGGLE TO DETERMINE SIMULATION/JOB NUMBER
# job = 1 # run the first job for example
# determine job number from pbs script
job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

################################################################################
# Set parameters that define the scenarios:

# random seed / simulation number
seed = tab$sim[tab$job == job]
# environment range of effect
env.range = tab$env_range[tab$job == job]
# latent range of effect
lat.range = tab$lat_range[tab$job == job]
# bias field range of effect
bias.range = tab$bias_range[tab$job == job]
# intercept for the PO data
int_po <- tab$int_po[tab$job == job]
# intercept for the PA data
int_pa <- tab$int_pa[tab$job == job]

# NOTE #
# int_po =
# -6.685 gets ~ 50 presences on average
# -5.3250 gets ~ 200 presences on average
# -3.9525 gets ~ 800 presences on average
#
# int_pa =
# -3.825 gets ~ 50 presence at survey sites on average
# -2.007 gets ~ 200 presence at survey sites on average
# 1.054 gets ~ 800 presence at survey sites on average

################################################################################

source("sim_PO_PA_data.R")
structured_data <- sim_occurrence_data(Intercept_po = int_po,
                                       Intercept_pa = int_pa,
                                       sites.sampled = 1000,
                                       rseed = seed,
                                       env.covariate.type = "random_field",
                                       presence.only.observer.bias.covariate.type = "random_field",
                                       presence.only.observer.bias.covariate.range = bias.range,
                                       env.covariate.range = env.range,
                                       latent.range = lat.range,
                                       latent.field = T,
                                       plotting = F
)
unstructured_data <- attr(structured_data, "presence-only")
quad <- attr(structured_data, "truth.grid")
c(n_pa = sum(structured_data$present), n_po = nrow(unstructured_data), ratio = sum(structured_data$present) / nrow(unstructured_data))

# add a spatial constraint to the survey data ##################################
################################################################################
# Interpolate some covariate at x, y locations #################################
interp.covar <- function(x.loc, y.loc, covar.name, data = quad){
  
  # turn the quadrature into a spatial pixels data frame
  sp.quad <- sp::SpatialPixelsDataFrame(points = quad[,c("x", "y")], data = quad[ , !colnames(quad) %in% c("x", "y", "quad.size")])
  
  # turn coordinates into SpatialPoints object:
  spp = sp::SpatialPoints(data.frame(x = x.loc,y = y.loc)) 
  
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- sp::over(spp, sp.quad[ , covar.name])
  v[is.na(v)] = 0 # NAs are a problem! Remove them
  return(v[,1])
}
################################################################################
survey <- data.frame(x = runif(1000, 50, 100), y = runif(1000, 1, 100))

# Interpolate the various fields #

# interpolate the environmental variable
survey$env <- interp.covar(x = survey$x, y = survey$y, covar.name = "env")
# calculate the linear predictor
survey$eta.fixed <- interp.covar(x = survey$x, y = survey$y, covar.name = "eta.fixed")
# calculate the true abundance rate
survey$mu <- interp.covar(x = survey$x, y = survey$y, covar.name = "mu")
# add in the latent field
survey$xi <- interp.covar(x = survey$x, y = survey$y, covar.name = "xi")
# calculate the presence probability
survey$pres.prob <-  1 - exp(-survey$mu)
# sample the presence/absence response
survey$present <- rbinom(nrow(survey), 1, survey$pres.prob)

# # Plot an example ###
# plot.res <- 500
# png(filename = paste0(getwd(), "/Results_spatially_constrained/Figures/our_sc_sim_setup.png"), res = plot.res, width = 6.2 * plot.res, height = 3.5 * plot.res)
# par(mfrow = c(1,2), mar = c(0, 0, 0, 1))
# # plot(vec2im(log(quad$lambda), quad$x, quad$y))
# image(x = sort(unique(quad$x)),
#       y = sort(unique(quad$y)),
#       z = scampr:::vec2mat(log(quad$lambda), quad$x, quad$y),
#       col = plasma(256), xlab = "", ylab = "", main = "", bty = 'n',
#       axes = F, asp = 1
# )
# title(main = expression(paste(PO, " Data and Intensity, ", lambda)), line = -0.5)
# points(unstructured_data[,c("x", "y")], col = "darkred")
# par(mar = c(0, 1, 0, 0))
# # plot(vec2im(log(quad$mu), quad$x, quad$y))
# image(x = sort(unique(quad$x)),
#       y = sort(unique(quad$y)),
#       z = scampr:::vec2mat(log(quad$mu), quad$x, quad$y),
#       col = plasma(256), xlab = "", ylab = "", main = "", bty = 'n',
#       axes = F, asp = 1
# )
# title(main = expression(paste(PA, " Data and Abundance Rate, ", mu)), line = -0.5)
# points(survey[,c("x", "y")], pch = c(4,1)[survey$present + 1], col = c("lightblue","darkblue")[survey$present + 1])
# par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
# dev.off()

# add a presence identifier to the quadrature
quad$present <- 0

# get the spatially constrained domain for prediction
quad_sc <- quad[quad$x < 50, ]

# created stacked unstructured_data and quad to be used as the data for scampr models
dat.scampr <- rbind(unstructured_data, quad)

# PA only model ##############################################################

# fit the base model without SRE
pa0 <- scampr(formula = present ~ env, data = structured_data, include.sre = F, sre.approx = "laplace", model.type = "PA")
# perform the basis opt.
pa <- basis.search.pa(pa0, domain.data = quad[,c("x","y")], max.basis.functions = 600, return.model = T)
# predict the mean abundance rate of the prediction points
pa.pred <- predict(pa, newdata = quad)
pa.pred_sc <- predict(pa, newdata = quad_sc)
pa.pred_po <- pa.pred

# fit the spatially constrained model
pa0_sc <- scampr(formula = present ~ env, data = survey, include.sre = F, sre.approx = "laplace", model.type = "PA")
pa_sc <- basis.search.pa(pa0_sc, domain.data = quad[,c("x","y")], max.basis.functions = 600, return.model = T)
pa_sc.pred <- predict(pa_sc, newdata = quad)
pa_sc.pred_sc <- predict(pa_sc, newdata = quad_sc)
pa_sc.pred_po <- pa_sc.pred

# PO only model ##############################################################

# fit the base model without SRE
po0 <- scampr(formula = present ~ env, data = dat.scampr, include.sre = F, sre.approx = "laplace", model.type = "PO")
# perform the basis opt.
po <- basis.search.po(po0, domain.data = quad[,c("x","y")], return.model = T, max.basis.functions = 600, which.approx = "laplace")
# predict the mean abundance rate of the prediction points
po.pred <- predict(po, newdata = quad)
po.pred_sc <- predict(po, newdata = quad_sc)
po.pred_po <- po.pred
  
# IDM ########################################################################

# fit the base model without SRE
idm0 <- scampr(present ~ env, data = dat.scampr, bias.formula = ~ 1, pa.data = structured_data, include.sre = F, sre.approx = "laplace", model.type = "IDM", latent.po.biasing = F)
idm <- basis.search(idm0, domain.data = quad[,c("x","y")], max.basis.functions = 600, return.model = T)
# predict the mean abundance rate of the prediction points
idm.pred <- predict(idm, newdata = quad)
idm.pred_sc <- predict(idm, newdata = quad_sc)
idm.pred_po <- predict(idm, newdata = quad, include.bias.accounting = T)

# fit the spatialy constrained model
idm0_sc <- scampr(present ~ env, data = dat.scampr, bias.formula = ~ 1, pa.data = survey, include.sre = F, sre.approx = "laplace", model.type = "IDM", latent.po.biasing = F)
idm_sc <- basis.search(idm0_sc, domain.data = quad[,c("x","y")], max.basis.functions = 600, return.model = T)
idm_sc.pred <- predict(idm_sc, newdata = quad)
idm_sc.pred_sc <- predict(idm_sc, newdata = quad_sc)
idm_sc.pred_po <- predict(idm_sc, newdata = quad, include.bias.accounting = T)

##############################################################################

# calculate metrics and return result
calc_KLdiv <- function(scampr.pred, domain.data) {
  m.prd <- exp(scampr.pred)
  return(as.numeric(domain.data$quad.size %*% (domain.data$mu * log(domain.data$mu / m.prd))) - as.numeric(domain.data$quad.size %*% (domain.data$mu - m.prd)))
}
calc_KLdiv_po <- function(scampr.pred, domain.data) {
  m.prd <- exp(scampr.pred)
  return(as.numeric(domain.data$quad.size %*% (domain.data$lambda * log(domain.data$lambda / m.prd))) - as.numeric(domain.data$quad.size %*% (domain.data$lambda - m.prd)))
}
calc_MAE <- function(scampr.pred, domain.data) {
  m.prd <- exp(scampr.pred)
  return(mean(abs((m.prd - mean(m.prd)) - (domain.data$mu - mean(domain.data$mu)))))
}

# set up the result data frame
res_po = data.frame(k = if(is.null(po$random.effects)){0}else{nrow(po$random.effects)}, k_bias = NA, Model = "PO Only", LL = logLik(po), KLdiv_full = calc_KLdiv(po.pred, quad), KLdiv_sc = NA, KLdiv_po = calc_KLdiv_po(po.pred_po, quad), MAE_full = calc_MAE(po.pred, quad), MAE_sc = NA, fit_flag = po$convergence, beta_est = po$fixed.effects["env", 1], beta_se = po$fixed.effects["env", 2])
res_pa = data.frame(k = if(is.null(pa$random.effects)){0}else{nrow(pa$random.effects)}, k_bias = NA, Model = "PA Only FULL", LL = logLik(pa), KLdiv_full = calc_KLdiv(pa.pred, quad), KLdiv_sc = calc_KLdiv(pa.pred_sc, quad_sc), KLdiv_po = calc_KLdiv_po(pa.pred_po, quad), MAE_full = calc_MAE(pa.pred, quad), MAE_sc = calc_MAE(pa.pred_sc, quad_sc), fit_flag = pa$convergence, beta_est = pa$fixed.effects["env", 1], beta_se = pa$fixed.effects["env", 2])
res_popa = data.frame(k = if(is.null(idm$random.effects)){0}else{nrow(idm$random.effects)}, k_bias = if(is.null(idm$random.bias.effects)){0}else{nrow(idm$random.bias.effects)}, Model = "IDM FULL", LL = logLik(idm), KLdiv_full = calc_KLdiv(idm.pred, quad), KLdiv_sc = calc_KLdiv(idm.pred_sc, quad_sc), KLdiv_po = calc_KLdiv_po(idm.pred_po, quad), MAE_full = calc_MAE(idm.pred, quad), MAE_sc = calc_MAE(idm.pred_sc, quad_sc), fit_flag = idm$convergence, beta_est = idm$fixed.effects["env", 1], beta_se = idm$fixed.effects["env", 2])
res_pa_sc = data.frame(k = if(is.null(pa_sc$random.effects)){0}else{nrow(pa_sc$random.effects)}, k_bias = NA, Model = "PA Only SC", LL = logLik(pa), KLdiv_full = calc_KLdiv(pa_sc.pred, quad), KLdiv_sc = calc_KLdiv(pa_sc.pred_sc, quad_sc), KLdiv_po = calc_KLdiv_po(pa_sc.pred_po, quad), MAE_full = calc_MAE(pa_sc.pred, quad), MAE_sc = calc_MAE(pa_sc.pred_sc, quad_sc), fit_flag = pa_sc$convergence, beta_est = pa_sc$fixed.effects["env", 1], beta_se = pa_sc$fixed.effects["env", 2])
res_popa_sc = data.frame(k = if(is.null(idm_sc$random.effects)){0}else{nrow(idm_sc$random.effects)}, k_bias = if(is.null(idm_sc$random.bias.effects)){0}else{nrow(idm_sc$random.bias.effects)}, Model = "IDM SC", LL = logLik(idm_sc), KLdiv_full = calc_KLdiv(idm_sc.pred, quad), KLdiv_sc = calc_KLdiv(idm_sc.pred_sc, quad_sc), KLdiv_po = calc_KLdiv_po(idm_sc.pred_po, quad), MAE_full = calc_MAE(idm_sc.pred, quad), MAE_sc = calc_MAE(idm_sc.pred_sc, quad_sc), fit_flag = idm_sc$convergence, beta_est = idm_sc$fixed.effects["env", 1], beta_se = idm_sc$fixed.effects["env", 2])

# collate the results
res_tab <- rbind(res_pa, res_po, res_popa, res_pa_sc, res_popa_sc)
res_tab$sim <- seed
res_tab$n_po <- nrow(unstructured_data)
res_tab$n_pa_pres <- sum(structured_data$present)
res_tab$bias_range = bias.range
res_tab$env_range = env.range
res_tab$latent_range = lat.range
res_tab$scen_pa = tab$scen_pa[tab$job == job]
res_tab$scen_po = tab$scen_po[tab$job == job]

# save the simulation result table in folder "Results/raw" inside the base dir
save(list = "res_tab", file = paste0(getwd(), "/Results_spatially_constrained/raw/res_", job, ".RDATA"))

