# INLA
home.wd <- getwd()

## Install required packages ###################################################
if(!require(scampr, quietly = T)){
  # scampr package can be installed from source provided in the code zip
  setwd("..")
  install.packages(paste0(getwd(), "/scampr_0.0.0.9000.tar.gz"), repos = NULL, type="source")
  library(scampr)
  setwd(home.wd)
}
if(!require(disdat, quietly = T)){
  install.packages("disdat")
  library(disdat)
}
if(!require(MASS, quietly = T)){
  install.packages("MASS")
  library(MASS)
}
if(!require(pROC, quietly = T)){
  install.packages("pROC")
  library(pROC)
}
if(!require(dplyr, quietly = T)){
  install.packages("dplyr")
  library(dplyr)
}
if(!require(ggplot2, quietly = T)){
  install.packages("ggplot2")
  library(ggplot2)
}
################################################################################

# Perform the spatial k-fold cross-validation #

# Here we are interested in demonstrating the benefit of the IDM for a single species
job = 7 # job = 2 doesn't select a random effects

# some hard-coded info
r = "NSW"
r_size <- 76.18 * 1000 # in meters^2
# will additionally include ordinal factors as categoricals in NSW
categoricalvars <- c("vegsys", "disturb", "soilfert")
flora_groups <- c("ot", "ou", "rt", "ru") # open-forest trees, open-forest understorey plants, rainforest trees, rainforest understorey plants

# reading presence-only and background species data for this region, one file per region:
presences <- disPo(r)
background <- disBg(r)
# reduce to just species of flora
presences <- presences[presences$group %in% flora_groups, ]

# get the relevant predictors
preds <- disPredictors(r)
preds <- preds[-13] # remove vegsys as this factor is likely correlated with the response (as it broadly categorises the main plant types in the area)

# names for all species - from supplementary material associated with Elith el al 2020
species <- unique(presences$spid)
# load("sample_sizes.RDATA")
# sp.numbers$job = match(sp.numbers$spid, species)

# obtain the presence/absence data
pa_list <- list()
for (grp in flora_groups) {
  pa_list[[grp]] <- cbind(disEnv(r, grp), disPa(r, grp)[-(1:4)])
}

# get the full domain data grid
load("NSW.RDATA") # load the pre-prepared full domain data
covars <- preds[!preds %in% c("disturb", "soilfert")]
covar_means <- colMeans(dat[ , covars], na.rm = T)
covar_sds <- sapply(dat[ , covars], sd, na.rm = T)
# calculate the scaled variables
tmp_covars_scaled <- scale(dat[, covars], center = covar_means, scale = covar_sds)
# retain the original variables
colnames(dat)[colnames(dat) %in% covars] <- paste0(colnames(dat)[colnames(dat) %in% covars], "_og")
# add in the scaled variables
dat <- cbind(dat, tmp_covars_scaled)
rm(tmp_covars_scaled)
gc()

## Model the particular species ##

# obtain the species
s = species[job]

# subset presence records of species for this species
sp_presence <- presences[presences$spid == s, ]
# add background data
dat_po <- rbind(sp_presence, background)
# add in the quadrature weights
dat_po$quad.size <- r_size / nrow(background)
dat_po$quad.size[dat_po$occ == 1] <- 0 # set quadrature weight to zero at the presence records

# identify the group of the species
grp <- sp_presence[, "group"][1]

# select out the corresponding PA dataset
dat_pa <- pa_list[[grp]]
# rename the species response variable to match the PO data (for scampr models)
dat_pa$occ <- dat_pa[ , s]
rm(presences, background, sp_presence, pa_list)
gc()

for(i in preds){
  if(i %in% categoricalvars){
    fac_col <- i
    if (fac_col == "soilfert") { # combining the soil fertility ratings beyond 3
      dat_po[ ,fac_col][dat_po[ ,fac_col] %in% c(4,5)] <- 3
      dat_pa[ ,fac_col][dat_pa[ ,fac_col] %in% c(4,5)] <- 3
      dat[ ,fac_col][dat[ ,fac_col] %in% c(4,5)] <- 3
    }
    dat[ ,fac_col] <- as.factor(dat[ ,fac_col])
    dat_po[ ,fac_col] <- as.factor(dat_po[ ,fac_col])
    dat_pa[ ,fac_col] <- as.factor(dat_pa[ ,fac_col])
    # expand and relevel incase the PO or PA datasets are missing levels from the full domain
    dat_po[ ,fac_col] <- forcats::fct_expand(dat_po[,fac_col], levels(dat[,fac_col]))
    dat_po[ ,fac_col] <- forcats::fct_relevel(dat_po[,fac_col], levels(dat[,fac_col]))
    dat_pa[ ,fac_col] <- forcats::fct_expand(dat_pa[,fac_col], levels(dat[,fac_col]))
    dat_pa[ ,fac_col] <- forcats::fct_relevel(dat_pa[,fac_col], levels(dat[,fac_col]))
  } else {
    # scale according to full domain means/sds
    dat_po[ , i] <- as.vector(scale(dat_po[ , i], center = covar_means[i], scale = covar_sds[i]))
    dat_pa[ , i] <- as.vector(scale(dat_pa[ , i], center = covar_means[i], scale = covar_sds[i]))
  }
}
gc()

# need to remove rows with NA predictor values
domain <- na.omit(dat)
rm(dat)
gc()

### Fit the models #############################################################

# set the model formula
form <- as.formula(paste0("occ ~ ", paste(preds, collapse = " + "), " + ", paste(paste0("I(", preds[!preds %in% c("disturb", "soilfert")], "^2)"), collapse = " + ")))

library(INLA)

# set the mesh according to the domain points
# mesh.time <- system.time(assign("mesh", inla.mesh.2d(loc.domain = domain[ , c("x","y")], max.edge=c(0.35,0.5), cutoff=0.05)))
mesh <- inla.mesh.2d(loc.domain = domain[ , c("x","y")], max.edge=c(0.35,0.5), cutoff=0.05)

# set the spde representation to be the mesh with defaults for parameters
# spde <- inla.spde2.matern(mesh)
spde <- inla.spde2.pcmatern(mesh,
                            prior.sigma = c(5, 0.01),
                            prior.range = c(0.1, 0.01)
)

# make A matrix for structured data
structured_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(dat_pa[ , c("x","y")]))

# make A matrix for unstructured data
unstructured_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(dat_po[dat_po$occ == 1 , c("x","y")]))

# make A matrix for quadrature points
quad_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(dat_po[dat_po$occ == 0 , c("x","y")]))

# make A matrix for the prediction points
pred_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(domain[ , c("x","y")]))

# Joint model

# One spatial field for shared effects and one for PO biasing
# Uses Simpson approach for PP data
# Binomial model for PA data
# Using cloglog

# set the number of integration points, presence points in PO data and prediction points
nq <- sum(dat_po$occ == 0)
n <- sum(dat_po$occ == 1)
np <- nrow(domain)


# change data to include 0s for nodes and 1s for presences
y.pp <- dat_po$occ

# add expectation vector (area for integration points/nodes and 0 for presences)
e.pp <- dat_po$quad.size

# make the full A matrix for the PPM component
A.pp <- rbind(unstructured_data_A, quad_data_A)

# need to make our own factor contrasts for INLA
fixed.po <- data.frame(get.design.matrix(form, dat_po))
fixed.pa <- data.frame(get.design.matrix(form, dat_pa))
fixed.pred <- data.frame(get.design.matrix(form, domain))

# unstructured data stack with integration points
stk_unstructured_data <- inla.stack(data=list(y=cbind(y.pp, NA), e = e.pp),
                                    effects=list(list(data.frame(interceptB=rep(1,nq+n)),
                                                      cti = dat_po$cti,
                                                      disturb2 = fixed.po$disturb2,
                                                      disturb3 = fixed.po$disturb3,
                                                      disturb4 = fixed.po$disturb4,
                                                      mi = dat_po$mi,
                                                      rainann = dat_po$rainann,
                                                      raindq = dat_po$raindq,
                                                      rugged = dat_po$rugged,
                                                      soildepth = dat_po$soildepth,
                                                      soilfert2 = fixed.po$soilfert2,
                                                      soilfert3 = fixed.po$soilfert3,
                                                      solrad = dat_po$solrad,
                                                      tempann = dat_po$tempann,
                                                      tempmin = dat_po$tempmin,
                                                      topo = dat_po$topo
                                    ),
                                    list(uns_field=1:spde$n.spde, bias_field = 1:spde$n.spde)),
                                    A=list(1,A.pp),
                                    tag="po_data")

# stack for structured data
# note intercept with different name
stk_structured_data <- inla.stack(data=list(y=cbind(NA, dat_pa$occ), Ntrials = rep(1, nrow(dat_pa))),
                                  effects=list(list(data.frame(interceptA=rep(1,nrow(dat_pa))),
                                                    cti = dat_pa$cti,
                                                    disturb2 = fixed.pa$disturb2,
                                                    disturb3 = fixed.pa$disturb3,
                                                    disturb4 = fixed.pa$disturb4,
                                                    mi = dat_pa$mi,
                                                    rainann = dat_pa$rainann,
                                                    raindq = dat_pa$raindq,
                                                    rugged = dat_pa$rugged,
                                                    soildepth = dat_pa$soildepth,
                                                    soilfert2 = fixed.pa$soilfert2,
                                                    soilfert3 = fixed.pa$soilfert3,
                                                    solrad = dat_pa$solrad,
                                                    tempann = dat_pa$tempann,
                                                    tempmin = dat_pa$tempmin,
                                                    topo = dat_pa$topo
                                  ),
                                  list(str_field=1:spde$n.spde)),
                                  A=list(1,structured_data_A),
                                  tag="pa_data")


# create the prediction stack
stk_pred_response <- inla.stack(data=list(y=cbind(rep(NA, np), rep(NA, np))),
                                effects = list(list(data.frame(interceptA=rep(1,np)),
                                                    cti = domain$cti,
                                                    disturb2 = fixed.pred$disturb2,
                                                    disturb3 = fixed.pred$disturb3,
                                                    disturb4 = fixed.pred$disturb4,
                                                    mi = domain$mi,
                                                    rainann = domain$rainann,
                                                    raindq = domain$raindq,
                                                    rugged = domain$rugged,
                                                    soildepth = domain$soildepth,
                                                    soilfert2 = fixed.pred$soilfert2,
                                                    soilfert3 = fixed.pred$soilfert3,
                                                    solrad = domain$solrad,
                                                    tempann = domain$tempann,
                                                    tempmin = domain$tempmin,
                                                    topo = domain$topo
                                ),
                                list(uns_field=1:spde$n.spde)),
                                A=list(1,pred_data_A),
                                tag='pred_response')

# combine the stacks
stk <- inla.stack(stk_unstructured_data, stk_structured_data, stk_pred_response)

# fit the model
result <- inla(y ~ interceptA + interceptB + cti + disturb2 + disturb3 + disturb4 + mi + rainann + raindq + rugged + soildepth + 
                 soilfert2 + soilfert3 + solrad + tempann + tempmin + topo + I(cti^2) + 
                 I(mi^2) + I(rainann^2) + I(raindq^2) + I(rugged^2) + I(soildepth^2) + 
                 I(solrad^2) + I(tempann^2) + I(tempmin^2) + I(topo^2) + 
                 f(uns_field, model = spde) + f(str_field, copy = "uns_field", fixed = TRUE) + f(bias_field, model = spde) -1,
               family=c("poisson", "binomial"),
               data=inla.stack.data(stk),
               control.predictor=list(A=inla.stack.A(stk), compute=TRUE),
               control.family = list(list(link = "log"),
                                     list(link = "cloglog")),
               E = inla.stack.data(stk)$e,
               Ntrials = inla.stack.data(stk)$Ntrials,
               control.compute = list(cpo=TRUE, waic = TRUE, dic = TRUE)
)

# create index to extract predictions
index.pred.response <- inla.stack.index(stk, tag="pred_response")$data
# add in the predictions
result$pred <- result$summary.fitted.values$mean[index.pred.response]

save(result, file = "INLA_MODEL.RDATA")
