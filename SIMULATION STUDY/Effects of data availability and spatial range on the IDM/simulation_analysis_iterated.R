if(!require(scampr, quietly = T)){
  # scampr package can be installed from github using devtools (see https://www.r-project.org/nosvn/pandoc/devtools.html for devtools installation)
  devtools::install_github("ElliotDovers/scampr", dependencies = "Imports", upgrade = "never")
  library(scampr)
}

# Get the job array
tab <- read.csv("job_array.csv")

# TOGGLE TO DETERMINE SIMULATION/JOB NUMBER
job = 1 # run the first job for example
# determine job number from pbs script
# job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

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

## for scampr models ###########################################################
library(scampr)

# add a presence identifier to the quadrature
quad$present <- 0

# created stacked unstructured_data and quad to be used as the data for scampr models
dat.scampr <- rbind(unstructured_data, quad)

# PA only model ##############################################################

# fit the base model without SRE
pa0 <- scampr(formula = present ~ env, data = structured_data, include.sre = F, sre.approx = "laplace", model.type = "PA")
# perform the basis opt.
pa <- basis.search.pa(pa0, domain.data = quad[,c("x","y")], max.basis.functions = 600, return.model = T)
# predict the mean abundance rate of the prediction points
pa_pred.time <- system.time(assign("pa.pred", predict(pa, newdata = quad)))
# collate the timings
pa.times <- as.numeric(pa$cpu["basis.search"] + pa_pred.time[3])

# PO only model ##############################################################

# fit the base model without SRE
po0 <- scampr(formula = present ~ env, data = dat.scampr, include.sre = F, sre.approx = "laplace", model.type = "PO")
# perform the basis opt.
po <- basis.search.po(po0, domain.data = quad[,c("x","y")], return.model = T, max.basis.functions = 600, which.approx = "laplace") # VA approx is fine here as we are going for speed!
# predict the mean abundance rate of the prediction points
po_pred.time <- system.time(assign("po.pred", predict(po, newdata = quad)))
# collate the timings
po.times <- as.numeric(po$cpu["basis.search"] + po_pred.time[3])

# IDM ########################################################################

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
idm_pred.time <- system.time(assign("idm.pred", predict(idm, newdata = quad)))
# collate the timings
idm.times <- as.numeric(pa$cpu["basis.search"] + po$cpu["basis.search"] + idm0$cpu["opt"] + idm$cpu["opt"] + idm_pred.time[3])

##############################################################################

# calculate metrics and return result
calc_KLdiv <- function(scampr.pred) {
  m.prd <- exp(scampr.pred)
  return(as.numeric(quad$quad.size %*% (quad$mu * log(quad$mu / m.prd))) - as.numeric(quad$quad.size %*% (quad$mu - m.prd)))
}
calc_MAE <- function(scampr.pred) {
  m.prd <- exp(scampr.pred)
  return(mean(abs((m.prd - mean(m.prd)) - (quad$mu - mean(quad$mu)))))
}

# set up the result data frame
res_po = data.frame(k = if(is.null(po$random.effects)){0}else{nrow(po$random.effects)}, k_bias = NA, Model = "PO Only", LL = logLik(po), KLdiv = calc_KLdiv(po.pred), MAE = calc_MAE(po.pred), timing = po.times, fit_flag = po$convergence)
res_pa = data.frame(k = if(is.null(pa$random.effects)){0}else{nrow(pa$random.effects)}, k_bias = NA, Model = "PA Only", LL = logLik(pa), KLdiv = calc_KLdiv(pa.pred), MAE = calc_MAE(pa.pred), timing = pa.times, fit_flag = pa$convergence)
res_popa = data.frame(k = if(is.null(idm$random.effects)){0}else{nrow(idm$random.effects)}, k_bias = if(is.null(idm$random.bias.effects)){0}else{nrow(idm$random.bias.effects)}, Model = "IDM", LL = logLik(idm), KLdiv = calc_KLdiv(idm.pred), MAE = calc_MAE(idm.pred), timing = idm.times, fit_flag = idm$convergence)

# collate the results
res_tab <- rbind(res_pa, res_po, res_popa)
res_tab$sim <- seed
res_tab$n_po <- nrow(unstructured_data)
res_tab$n_pa_pres <- sum(structured_data$present)
res_tab$bias_range = bias.range
res_tab$env_range = env.range
res_tab$latent_range = lat.range
res_tab$scen_pa = tab$scen_pa[tab$job == job]
res_tab$scen_po = tab$scen_po[tab$job == job]

# save the simulation result table in folder "Results/raw" inside the base dir
save(list = "res_tab", file = paste0(getwd(), "/Results/raw/res_", job, ".RDATA"))

