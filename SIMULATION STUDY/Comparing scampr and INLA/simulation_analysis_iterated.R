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
if(!require(rgeos, quietly = T)){
  install.packages("rgeos")
  library(rgeos)
}
if(!require(deldir, quietly = T)){
  install.packages("deldir")
  library(deldir)
}
################################################################################

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
env.range = 15
# latent range of effect
lat.range = 25
# bias field range of effect
bias.range = 29
# model
model_to_test <- tab$fit_model[tab$job == job]
# quadrature
quad_to_use <- tab$quad[tab$job == job]
# number of predictions
n_pred_to_use <- tab$n_pred[tab$job == job]

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
# c(n_pa = sum(structured_data$present), n_po = nrow(unstructured_data), ratio = sum(structured_data$present) / nrow(unstructured_data))

################################################################################
# create the INLA mesh to be used. Max edge length of 10 units should be fine enough to model the latent field (25 units) and bias field (20 units)
mesh <- inla.mesh.2d(loc.domain = domain.data[ , c("x", "y")], max.edge=c(10,30), cutoff=2, offset = c(5,20))
# get the boundary for the domain
bnd <- t(matrix(c(min(domain.data$x),min(domain.data$x),max(domain.data$x),min(domain.data$x),max(domain.data$x),max(domain.data$y),min(domain.data$y),max(domain.data$y),min(domain.data$y),min(domain.data$y)), 2))
# make dual mesh
dd <- deldir::deldir(mesh$loc[, 1], mesh$loc[, 2])
tiles <- deldir::tile.list(dd)
# make domain into spatial polygon
domainSP <- SpatialPolygons(list(Polygons(
  list(Polygon(bnd)), '0')))
# intersection between domain and dual mesh
poly.gpc <- as(domainSP@polygons[[1]]@Polygons[[1]]@coords, "gpc.poly")
# w now contains area of voronoi polygons
w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x, p$y), "gpc.poly"), poly.gpc)))
################################################################################

# set the quadrature according to the job (full domain data or the INLA mesh, n_quad = 10,000 or 415)
if (quad_to_use == "full") {
  quad <- domain.data
  quad.in <- rep(T, nrow(quad))
} else {
  quad <- data.frame(mesh$loc[ , 1:2])
  colnames(quad) <- c("x", "y")
  # find quad point inside the boundary
  quad.in <- in.out(bnd, as.matrix(quad))
  # set the quadrature sizes
  quad$quad.size <- w
  # interpolate all the continuous variables
  for (colm in colnames(domain.data)[!colnames(domain.data) %in% c("x", "y", "quad.size")]) {
    eval(parse(text = paste0("quad$", colm, " <- interp.covar(quad$x, quad$y, '", colm, "', domain.data)")))
  }
  # # visual sense check
  # colfn <- function(x, colors=plasma(100), colsteps=100) {
  #   return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
  # }
  # par(mfrow = c(2,1))
  # plot(tiles, fillcol = colfn(quad$env), bty = "n", main = "Reduced Resolution")
  # plot(vec2im(domain.data$env, domain.data$x, domain.data$y), main = "Fine scale Truth")
  # par(mfrow = c(1,1))
}

# set the prediction locations according to the job (full domain data or the INLA mesh, n_pred = 10,000 or 100)
if (n_pred_to_use == "full") {
  pred <- domain.data
  pred.in <- rep(T, nrow(quad))
} else {
  resolution <- c(10,10)
  pred <- expand.grid(x=seq(resolution[1]/2, max(quad$x), resolution[1]), 
                      y=seq(resolution[2]/2, max(quad$y), resolution[2]))
  # get the quadrature sizes
  dx <- unique(diff(sort(unique(pred$x))))[1] # [1] required due to machine error
  dy <- unique(diff(sort(unique(pred$y))))[1] # [1] required due to machine error
  D.size <- diff(range(pred$x)) * diff(range(pred$y))
  ## Quadrature Sizes ##
  # rectangles centered at the quadrature
  pred$quad.size <- rep(dx * dy, nrow(pred))
  # along the edges of the lower square we are missing a unit in one direction
  pred$quad.size[pred$x == min(pred$x) | pred$y == min(pred$y)] <- (dx - 1) * dy
  # small lower left corner is missing one unit from both
  pred$quad.size[pred$x == min(pred$x) & pred$y == min(pred$y)] <- (dx - 1) * (dy - 1)
  # interpolate all the continuous variables
  for (colm in colnames(domain.data)[!colnames(domain.data) %in% c("x", "y", "quad.size")]) {
    eval(parse(text = paste0("pred$", colm, " <- interp.covar(pred$x, pred$y, '", colm, "', domain.data)")))
  }
  # interpolation can introduce zero mean abundances which we will set to the minimum within the domain data
  pred$mu[pred$mu <= 0] <- min(domain.data$mu)
  
  ### IF WE WANT TO PREDICT AT THE MESH VERTICES
  # pred <- data.frame(mesh$loc[ , 1:2])
  # colnames(pred) <- c("x", "y")
  # # find prediction points inside the boundary
  # pred.in <- in.out(bnd, as.matrix(pred))
  # # set the quadrature sizes
  # pred$quad.size <- w
  # # interpolate all the continuous variables
  # for (colm in colnames(domain.data)[!colnames(domain.data) %in% c("x", "y", "quad.size")]) {
  #   eval(parse(text = paste0("pred$", colm, " <- interp.covar(pred$x, pred$y, '", colm, "', domain.data)")))
  # }
  # # reduce the prediction locations to within the domain
  # pred <- pred[pred.in, ]
  # # interpolation can introduce zero mean abundances which we will set to the minimum within the domain data
  # pred$mu[pred$mu <= 0] <- min(domain.data$mu)
}

# fit and predict using the appropriate model for the job
if (model_to_test == "INLA PA") {
  source("inla_pa.R")
  res <- inla_pa(structured_data = structured_data, mesh = mesh, pred = pred)
} else if (model_to_test == "INLA PO") {
  source("inla_po.R")
  res <- inla_po(unstructured_data = unstructured_data, quad = quad, mesh = mesh, pred = pred)
} else if (model_to_test == "INLA IDM") {
  source("inla_idm.R")
  res <- inla_idm(structured_data = structured_data, unstructured_data = unstructured_data, quad = quad, mesh = mesh, pred = pred)
} else if (model_to_test == "SCAMPR") {
  source("scampr_all.R")
  res <- scampr_all(structured_data = structured_data, unstructured_data = unstructured_data, quad = quad, pred = pred, domain.data = domain.data, prune.n = 4)
} else {
  source("scampr_fixed_all.R")
  res <- scampr_fixed_all(structured_data = structured_data, unstructured_data = unstructured_data, quad = quad, pred = pred, domain.data = domain.data, prune.n = 4)
}

# collate with sim info
res_tab <- cbind(tab[tab$job == job, ], res)
if (model_to_test == "SCAMPR") {
  attr(res_tab, "all") <- attr(res, "all")
}
# adjust the parameters
res_tab$env_range <- env.range
res_tab$bias_range <- bias.range
res_tab$lat_range <- lat.range
res_tab$n_po <- nrow(unstructured_data)
res_tab$n_pres_pa <- sum(structured_data$present)

# save the simulation result table in folder "Results" inside the base dir
save(list = "res_tab", file = paste0(getwd(), "/Results/res_", job, ".RDATA"))