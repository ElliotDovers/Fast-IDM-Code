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
if(!require(INLA, quietly = T)){
  install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  library(INLA)
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
if(!require(mgcv, quietly = T)){
  install.packages("mgcv")
  library(mgcv)
}
# if(!require(rgeos, quietly = T)){
#   install.packages("rgeos")
#   library(rgeos)
# }
# if(!require(deldir, quietly = T)){
#   install.packages("deldir")
#   library(deldir)
# }
################################################################################

# Get the job array
tab <- read.csv("job_array_all.csv")

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
# model
model_to_test <- tab$fit_model[tab$job == job]

##############################################################################
# Interpolate some covariate at x, y locations ###############################
interp.covar <- function(x.loc, y.loc, covar.name, domain.data){
  
  # turn the quadrature into a spatial pixels data frame
  sp.domain <- sp::SpatialPixelsDataFrame(points = domain.data[,c("x", "y")], data = domain.data[ , !colnames(domain.data) %in% c("x", "y", "quad.size")])
  
  # turn coordinates into SpatialPoints object:
  spp = sp::SpatialPoints(data.frame(x = x.loc,y = y.loc)) 
  
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- sp::over(spp, sp.domain[ , covar.name])
  v[is.na(v)] = 0 # NAs are a problem! Remove them
  return(v[,1])
}
################################################################################

################################################################################

# simulate the data
source("sim_PO_PA_data.R")
structured_data <- sim_occurrence_data(Intercept_po = -3.9525, # gets ~ 800 on average
                                       Intercept_pa = -2.007, # gets ~ 200 presence at survey sites
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
domain.data <- attr(structured_data, "truth.grid")

################################################################################
# create the INLA mesh to be used. Max edge length of 10 units should be fine enough to model the latent field (25 units) and bias field (20 units)
tmp.time <- system.time(assign("mesh", inla.mesh.2d(loc.domain = domain.data[ , c("x", "y")], max.edge=c(10,30), cutoff=2, offset = c(5,20))))
mesh$timing.init <- tmp.time[3]

# set the quadrature
quad <- domain.data

# set the prediction locations
pred <- domain.data

# fit and predict using the appropriate model for the job
if (model_to_test == "INLA") {
  source("inla_all.R")
  res <- inla_all(structured_data = structured_data, unstructured_data = unstructured_data, quad = quad, mesh = mesh, pred = pred)
} else if (model_to_test == "SCAMPR") {
  source("scampr_all.R")
  res <- scampr_all(structured_data = structured_data, unstructured_data = unstructured_data, quad = quad, pred = pred, domain.data = domain.data, prune.n = 4)
} else if (model_to_test == "SCAMPR FIXED") {
  source("scampr_fixed_all.R")
  res <- scampr_fixed_all(structured_data = structured_data, unstructured_data = unstructured_data, quad = quad, pred = pred, domain.data = domain.data, prune.n = 4)
} else if (model_to_test == "MGCV") {
  source("mgcv_all.R")
  res <- mgcv_all(structured_data = structured_data, unstructured_data = unstructured_data, quad = quad, pred = pred, domain.data = domain.data)
} else {
  source("mgcv_fixed_all.R")
  res <- mgcv_fixed_all(structured_data = structured_data, unstructured_data = unstructured_data, quad = quad, pred = pred, domain.data = domain.data)
}

# collate with sim info
res_tab <- cbind(tab[tab$job == job, ], res)

# adjust the parameters
res_tab$env_range <- env.range
res_tab$bias_range <- bias.range
res_tab$lat_range <- lat.range
res_tab$n_po <- nrow(unstructured_data)
res_tab$n_pres_pa <- sum(structured_data$present)

# save the simulation result table in folder "Results" inside the base dir
save(list = "res_tab", file = paste0(getwd(), "/Results_all/res_", job, ".RDATA"))