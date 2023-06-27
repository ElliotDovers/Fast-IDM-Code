home.wd <- getwd()

## Install required packages ###################################################
if(!require(scampr, quietly = T)){
  # scampr package can be installed from source provided in the code zip
  setwd("..")
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
if(!require(sp, quietly = T)){
  install.packages("sp")
  library(sp)
}
if(!require(fields, quietly = T)){
  install.packages("fields")
  library(fields)
}
if(!require(resahpe2, quietly = T)){
  install.packages("reshape2")
  library(reshape2)
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

# seed
seed = tab$Sim[tab$job == job]
# nsamp
structured_sample_size <- tab$N_pa[tab$job == job]
# lambda (intercept of the linear predictor)
lambda = -2
# env.beta (coefficient on the environmental covariate)
env.beta = 2
# dim (finest resolution of the truth grid)
dim = c(300,300)
# sigma2x
sigma2x = 0.5
# kappa
kappa = 0.05
# strata
strata = 25
# rows
rows = 5
# cols
cols = 5
# probs
probs = tab$Prob_PO[tab$job == job]
# qsize
qsize = 1
# rho
rho = 0.99
# resolution (of the reduced prediction grid)
resolution = c(10,10)
# correlation (between the env covariate and bias covariate)
correlation = tab$Correlated[tab$job == job] == 1
# quadrature points type
qtype <- "full" # set to "mesh" to use the INLA mesh as quadrature points

################################################################################

# removing so they can all have same truth
source("Functions to generate data and sample.R")
g1 <- genDataFunctions(dim = dim,
                       lambda = lambda,
                       env.beta = env.beta,
                       seed = seed,
                       kappa = kappa,
                       sigma2x = sigma2x,
                       strata = strata,
                       rows = rows,
                       cols = cols,
                       probs = probs,
                       nsamp = structured_sample_size,
                       plot = FALSE,
                       plotdat = FALSE,
                       qsize = qsize, 
                       rho = rho,
                       correlated = correlation)

structured_data <- g1$structured_data
unstructured_data <- g1$unstructured_data
biasfield <- g1$biasfield
dat1 <- g1$dat1
biascov <- g1$biascov
num_presence <- sum(structured_data$presence)
# NOTES:
# strata1$sim1 == biasfield$sim1... AND this is equal to an axis flipped version of dat1$rf.s!!!
# cbind(strata1$sim1[order(strata1$x)], biasfield$sim1[order(biasfield$x)], quad$abund) ARE ALL THE SAME (where quad$abund is vectorised dat1$rf.s)
# dat1$rf.s == dat1$Lam
# biascov == biasfield$covariate
# dat1$rf.s appears to be the true rate of occurrence which is then thinned by biasfield$stratprobs
# so maybe PO intensity is quad$abund + log(-log(1-biasfield$stratprobs))?
# since the full presences are thinned by rbinom(nrow(dat1$xy),1,biasfield$stratprobs)

## for scampr models ###########################################################
cov.vec <- NULL
bias.vec <- NULL
log_mu_A.vec <- NULL
tmp <- NULL
xs <- NULL
ys <- NULL
for (i in 1:dim[1]) {
  cov.vec[(1:dim[2]) + (i-1)*dim[2]] <- dat1$gridcov[i, ]
  bias.vec[(1:dim[2]) + (i-1)*dim[2]] <- biascov[i, ]
  log_mu_A.vec[(1:dim[2]) + (i-1)*dim[2]] <- dat1$rf.s[i, ] # rf.s == log(Lam) within dat1
  tmp[(1:dim[2]) + (i-1)*dim[2]] <- biascov[i, ]
  xs <- c(xs, 1:dim[1])
  ys <- c(ys, rep((1:dim[2])[i], dim[2]))
}
# the true presence intensity should be equal to the inverted biasfield sim1
all(biasfield$sim1[order(biasfield$x)] == log_mu_A.vec)
# calculate the thinned intensity that generates the presence-only data
log_lambda_po <- log_mu_A.vec + log(-log(1 - biasfield$stratprobs))

# set up the full quadrature
quad <- data.frame(log_mu_A=log_mu_A.vec,log_lambda_po=log_lambda_po,x=xs,y=ys,env=cov.vec,bias=bias.vec,presence=0,quad.size = 1)

# add in the true intensities to the unstructured data
unstructured_data$log_mu_A <- quad$log_mu_A[match(paste(unstructured_data$x, unstructured_data$y, sep = ":"), paste(quad$x, quad$y, sep = ":"))]
unstructured_data$log_lambda_po <- quad$log_lambda_po[match(paste(unstructured_data$x, unstructured_data$y, sep = ":"), paste(quad$x, quad$y, sep = ":"))]
# add in the quadrat sizes
unstructured_data$quad.size <- 0
# add in the true intensity to the structured data
structured_data$log_mu_A <- quad$log_mu_A[match(paste(structured_data$x, structured_data$y, sep = ":"), paste(quad$x, quad$y, sep = ":"))]

# adjust the quadrature type used for scampr
if (qtype == "mesh") {
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
  library(INLA)
  library(reshape2)
  library(deldir)
  library(rgeos)
  library(fields)
  mesh <- inla.mesh.2d(loc.domain = biasfield[,c(1,2)], max.edge=c(20,40), cutoff=2, offset = c(5,20))
  max_x <- max(biasfield$x)
  max_y <- max(biasfield$y)
  loc.d <- t(matrix(c(0,0,max_x,0,max_x,max_y,0,max_y,0,0), 2))
  dd <- deldir::deldir(mesh$loc[, 1], mesh$loc[, 2])
  tiles <- deldir::tile.list(dd)
  domainSP <- SpatialPolygons(list(Polygons(list(Polygon(loc.d)), '0')))
  poly.gpc <- as(domainSP@polygons[[1]]@Polygons[[1]]@coords, "gpc.poly")
  # w now contains area of voronoi polygons
  w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x,p$y), "gpc.poly"), poly.gpc)))
  
  # set up the quadrature for scampr
  s.quad <- data.frame(mesh$loc[w > 0 , 1:2]) # cut to those points in the domain (with positive weight)
  colnames(s.quad) <- c("x", "y")
  # set the quadrature sizes
  s.quad$quad.size <- w[w > 0]
  # interpolate all the continuous variables
  for (colm in colnames(quad)[!colnames(quad) %in% c("x", "y", "quad.size")]) {
    eval(parse(text = paste0("s.quad$", colm, " <- interp.covar(s.quad$x, s.quad$y, '", colm, "', domain.data = quad)")))
  }
} else {
  # set.seed(1)
  M <- 15000 # fix the number of quadrature points to sample
  quad.sample.id <- sample(1:nrow(quad), M, replace = F)
  # set.seed(NULL)
  s.quad <- quad[quad.sample.id, ]
  # add in quadrat sizes
  s.quad$quad.size <- prod(dim) / M
}

# created stacked unstructured_data and quad to be used as the data for scampr models
dat.scampr <- rbind(unstructured_data[ , c("x", "y", "quad.size", "presence", "env", "bias", "log_mu_A", "log_lambda_po")], s.quad[, c("x", "y", "quad.size", "presence", "env", "bias", "log_mu_A", "log_lambda_po")])

# set the basis functions
bfs <- scampr::simple_basis(20, data = quad) # k = 400
bfs_bias <- scampr::simple_basis(5, data = quad) # k = 25; need a second set of broad basis functions for the scampr version of 'jointtwo'
################################################################################

source("validation_function.R")

# Model A (PA data only)
source("Run models structured.R")
try(assign("m.a", structured_model(structured_data, dat1, biasfield, plotting = FALSE)))
if (exists("m.a")) {
  v.a <- validation_function(result=m.a[[2]], resolution=c(10,10), join.stack=m.a[[1]], mesh = m.a$mesh, model_type="structured", unstructured_data=unstructured_data,
                             structured_data = structured_data, dat1 = dat1, summary_results=T, qsize = 1, absolute=FALSE, plotting = FALSE, dim = dim)
} else {
  tab <- data.frame(Model = "structured", MAE = NA)
  attr(tab, "alt") <- data.frame(Model = "structured", KL_div = NA, KL_presprob = NA, KL_quad_po = NA, timing = NA, fit_flag = 1)
  v.a = list(tab, correlation = NA, coefficients = NA)
  names(v.a) <- c("Proto-table", "correlation", "coefficients")
  rm(tab)
}
m.a.s <- scampr::scampr(formula = presence ~ env, data = structured_data, basis.functions = bfs, sre.approx = "laplace", model.type = "PA")
v.a.s <- validation_function(result=m.a.s, resolution=c(10,10), model_type="structured", unstructured_data=unstructured_data, structured_data = structured_data, dat1 = dat1,
                             summary_results=T, qsize = 1, absolute=FALSE, plotting = FALSE, dim = dim, is.scampr = T)

# Model B (PO data only)
source("Run models.R")
try(assign("m.b", unstructured_model(unstructured_data, dat1, biasfield, dim = dim, plotting = FALSE)))
if (exists("m.b")) {
  v.b <- validation_function(result=m.b[[2]], resolution=c(10,10), join.stack=m.b[[1]], mesh = m.b$mesh, model_type="unstructured",
                             unstructured_data = unstructured_data, structured_data = structured_data, dat1 = dat1, summary_results=T, absolute=FALSE, dim = dim)
} else {
  tab <- data.frame(Model = "unstructured", MAE = NA)
  attr(tab, "alt") <- data.frame(Model = "unstructured", KL_div = NA, KL_presprob = NA, KL_quad_po = NA, timing = NA, fit_flag = 1)
  v.b = list(tab, correlation = NA, coefficients = NA)
  names(v.b) <- c("Proto-table", "correlation", "coefficients")
  rm(tab)
}
m.b.s <- scampr::scampr(presence ~ env, data = dat.scampr, basis.functions = bfs, sre.approx = "laplace", model.type = "PO")
v.b.s <- validation_function(result=m.b.s, resolution=c(10,10), model_type="unstructured", unstructured_data=unstructured_data,
                                    structured_data = structured_data, dat1 = dat1, summary_results=T, qsize = 1, absolute=FALSE, dim = dim, is.scampr = T)

# Model C (IDM without accounting for bias)
source("Run models joint.R")
try(assign("m.c", joint_model(structured_data, unstructured_data, dat1, biasfield, plotting = FALSE)))
if (exists("m.c")) {
  v.c <- validation_function(result=m.c[[2]], resolution=c(10,10), join.stack=m.c[[1]], mesh = m.c$mesh, model_type="joint",
                             unstructured_data = unstructured_data, structured_data = structured_data,
                             dat1 = dat1, summary_results=T, absolute=FALSE, dim = dim)
} else {
  tab <- data.frame(Model = "joint", MAE = NA)
  attr(tab, "alt") <- data.frame(Model = "joint", KL_div = NA, KL_presprob = NA, KL_quad_po = NA, timing = NA, fit_flag = 1)
  v.c = list(tab, correlation = NA, coefficients = NA)
  names(v.c) <- c("Proto-table", "correlation", "coefficients")
  rm(tab)
}
m.c.s <- scampr::scampr(presence ~ env, data = dat.scampr, bias.formula = ~ 1, IDM.presence.absence.df = structured_data, basis.functions = bfs, sre.approx = "laplace", model.type = "IDM")
v.c.s <- validation_function(result=m.c.s, resolution=c(10,10), model_type="joint", unstructured_data=unstructured_data,
                    structured_data = structured_data, dat1 = dat1, summary_results=T, qsize = 1, absolute=FALSE, dim = dim, is.scampr = T)

# Model D (PO data only with bias covariate)
source("Run models unstructured bias covariate.R")
try(assign("m.d", unstructured_model_cov(unstructured_data, dat1, biasfield, dim = dim, plotting = FALSE, biascov=biascov)))
if (exists("m.d")) {
  v.d <- validation_function(result=m.d[[2]], resolution=c(10,10), join.stack=m.d[[1]], mesh = m.d$mesh, model_type="unstructuredcov",
                             unstructured_data = unstructured_data, structured_data = structured_data, dat1 = dat1, summary_results=T, absolute=FALSE, dim = dim)
} else {
  tab <- data.frame(Model = "unstructuredcov", MAE = NA)
  attr(tab, "alt") <- data.frame(Model = "unstructuredcov", KL_div = NA, KL_presprob = NA, KL_quad_po = NA, timing = NA, fit_flag = 1)
  v.d = list(tab, correlation = NA, coefficients = NA)
  names(v.d) <- c("Proto-table", "correlation", "coefficients")
  rm(tab)
}
m.d.s <- scampr::scampr(presence ~ env, data = dat.scampr, bias.formula = ~ bias - 1, basis.functions = bfs, sre.approx = "laplace", model.type = "PO")
v.d.s <- validation_function(result=m.d.s, resolution=c(10,10), model_type="unstructuredcov", unstructured_data=unstructured_data,
                                    structured_data = structured_data, dat1 = dat1, summary_results=T, absolute=FALSE, dim = dim, is.scampr = T)

# Model E (IDM with bias covariate)
source("Run models joint covariate for bias.R")
try(assign("m.e", joint_model_cov(structured_data, unstructured_data, dat1, biasfield, biascov=biascov)))
if (exists("m.e")) {
  v.e <- validation_function(result=m.e[[2]], resolution=c(10,10), join.stack=m.e[[1]], mesh = m.e$mesh, model_type="jointcov",
                             unstructured_data = unstructured_data, structured_data = structured_data,
                             dat1 = dat1, summary_results=T, absolute = FALSE, dim = dim)
} else {
  tab <- data.frame(Model = "jointcov", MAE = NA)
  attr(tab, "alt") <- data.frame(Model = "jointcov", KL_div = NA, KL_presprob = NA, KL_quad_po = NA, timing = NA, fit_flag = 1)
  v.e = list(tab, correlation = NA, coefficients = NA)
  names(v.e) <- c("Proto-table", "correlation", "coefficients")
  rm(tab)
}
m.e.s <- scampr::scampr(presence ~ env, data = dat.scampr, bias.formula = ~ bias, IDM.presence.absence.df = structured_data, basis.functions = bfs, sre.approx = "laplace", model.type = "IDM")
v.e.s <- validation_function(result=m.e.s, resolution=c(10,10), model_type="jointcov", unstructured_data=unstructured_data,
                    structured_data = structured_data, dat1 = dat1, summary_results=T, absolute=FALSE, dim = dim, is.scampr = T)

# Model F (IDM with secondary latent field and without bias covariate)
source("Run models joint second spatial field.R")
try(assign("m.f", joint_model2(structured_data, unstructured_data, dat1, biasfield)))
if (exists("m.f")) {
  v.f <- validation_function(result=m.f[[2]], resolution=c(10,10), join.stack=m.f[[1]], mesh = m.f$mesh, model_type="jointtwo",
                             unstructured_data = unstructured_data, structured_data = structured_data,
                             dat1 = dat1, summary_results=T, absolute = FALSE, dim = dim)
} else {
  tab <- data.frame(Model = "jointtwo", MAE = NA)
  attr(tab, "alt") <- data.frame(Model = "jointtwo", KL_div = NA, KL_presprob = NA, KL_quad_po = NA, timing = NA, fit_flag = 1)
  v.f = list(tab, correlation = NA, coefficients = NA)
  names(v.f) <- c("Proto-table", "correlation", "coefficients")
  rm(tab)
}
m.f.s1 <- scampr::scampr(presence ~ env, data = dat.scampr, bias.formula = ~ 1, IDM.presence.absence.df = structured_data, basis.functions = bfs, sre.approx = "laplace", model.type = "IDM", latent.po.biasing = T, po.biasing.basis.functions = bfs_bias)
v.f.s1 <- validation_function(result=m.f.s1, resolution=c(10,10), model_type="jointtwo", unstructured_data=unstructured_data,
                             structured_data = structured_data, dat1 = dat1, summary_results=T, absolute=FALSE, dim = dim, is.scampr = T)
m.f.s2 <- scampr::scampr(presence ~ env, data = dat.scampr, bias.formula = ~ 1, IDM.presence.absence.df = structured_data, basis.functions = bfs_bias, sre.approx = "laplace", model.type = "IDM", latent.po.biasing = T, po.biasing.basis.functions = bfs)
v.f.s2 <- validation_function(result=m.f.s2, resolution=c(10,10), model_type="jointtwo", unstructured_data=unstructured_data,
                              structured_data = structured_data, dat1 = dat1, summary_results=T, absolute=FALSE, dim = dim, is.scampr = T)

# decide on the best configuration:
choice <- which.max(c(logLik(m.f.s1), logLik(m.f.s2)))
if (choice == 1) {
  m.f.s <- m.f.s1
  v.f.s <- v.f.s1
} else {
  m.f.s <- m.f.s2
  v.f.s <- v.f.s2
}
# update the timing to include that of both models
attr(v.f.s$`Proto-table`, "alt")$timing <- attr(v.f.s1$`Proto-table`, "alt")$timing + attr(v.f.s2$`Proto-table`, "alt")$timing

# Create data frame for saving individual sim results
res_tab <- data.frame(rbind(
cbind(Fit = "INLA", v.a$`Proto-table`, Correlation = v.a$correlation, attr(v.a$`Proto-table`, "alt")[-1], ll = m.a$result$mlik[2], search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "scampr", v.a.s$`Proto-table`, Correlation = v.a.s$correlation, attr(v.a.s$`Proto-table`, "alt")[-1], ll = logLik(m.a.s), search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "INLA", v.b$`Proto-table`, Correlation = v.b$correlation, attr(v.b$`Proto-table`, "alt")[-1], ll = m.b$result$mlik[2], search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "scampr", v.b.s$`Proto-table`, Correlation = v.b.s$correlation, attr(v.b.s$`Proto-table`, "alt")[-1], ll = logLik(m.b.s), search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "INLA", v.c$`Proto-table`, Correlation = v.c$correlation, attr(v.c$`Proto-table`, "alt")[-1], ll = m.c$result$mlik[2], search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "scampr", v.c.s$`Proto-table`, Correlation = v.c.s$correlation, attr(v.c.s$`Proto-table`, "alt")[-1], ll = logLik(m.c.s), search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "INLA", v.d$`Proto-table`, Correlation = v.d$correlation, attr(v.d$`Proto-table`, "alt")[-1], ll = m.d$result$mlik[2], search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "scampr", v.d.s$`Proto-table`, Correlation = v.d.s$correlation, attr(v.d.s$`Proto-table`, "alt")[-1], ll = logLik(m.d.s), search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "INLA", v.e$`Proto-table`, Correlation = v.e$correlation, attr(v.e$`Proto-table`, "alt")[-1], ll = m.e$result$mlik[2], search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "scampr", v.e.s$`Proto-table`, Correlation = v.e.s$correlation, attr(v.e.s$`Proto-table`, "alt")[-1], ll = logLik(m.e.s), search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "INLA", v.f$`Proto-table`, Correlation = v.f$correlation, attr(v.f$`Proto-table`, "alt")[-1], ll = m.f$result$mlik[2], search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "scampr", v.f.s$`Proto-table`, Correlation = v.f.s$correlation, attr(v.f.s$`Proto-table`, "alt")[-1], ll = logLik(m.f.s), search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "scampr_ll_k_dense_k_bias_coarse", v.f.s1$`Proto-table`, Correlation = v.f.s1$correlation, attr(v.f.s1$`Proto-table`, "alt")[-1], ll = logLik(m.f.s1), search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation),
cbind(Fit = "scampr_ll_k_coarse_k_bias_dense", v.f.s2$`Proto-table`, Correlation = v.f.s2$correlation, attr(v.f.s2$`Proto-table`, "alt")[-1], ll = logLik(m.f.s2), search.time = NA, sim = seed, n_po = nrow(unstructured_data), n_pa = structured_sample_size, detect_prob_po = probs, correlated_bias = correlation)
))

# save the simulation result table in folder "Results/raw" inside the base dir
save(list = "res_tab", file = paste0(getwd(), "/Results/raw/res_", job, ".RDATA"))

