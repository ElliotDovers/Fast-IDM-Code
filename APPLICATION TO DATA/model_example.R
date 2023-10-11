home.wd <- getwd()

# get functions to fit IDMs using scampr and INLA and make predictions across the domain
source("idm_scampr_w_prediction.R")
source("idm_inla_w_prediction.R")

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
load("sample_sizes.RDATA")
sp.numbers$job = match(sp.numbers$spid, species)

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

### Fit the models #############################################################

# set the model formula
form <- as.formula(paste0("occ ~ ", paste(preds, collapse = " + "), " + ", paste(paste0("I(", preds[!preds %in% c("disturb", "soilfert")], "^2)"), collapse = " + ")))

# fit the base scampr models
base_po <- scampr(form, data = dat_po, include.sre = F, model.type = "PO", sre.approx = "laplace")
base_pa <- scampr(form, data = dat_pa, include.sre = F, model.type = "PA", sre.approx = "laplace")
base_idm <- scampr(form, data = dat_po, bias.formula = ~ 1, IDM.presence.absence.df = dat_pa, include.sre = F, model.type = "IDM", sre.approx = "laplace", latent.po.biasing = F)

# fit the PA model
res_pa <- basis.search.pa(base_pa, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10, max.basis.functions = 200)
# fit the PO model
res_po <- basis.search.po(base_po, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10, which.approx = "laplace", max.basis.functions = 100)

# fit the IDM - with the function that also does predictions (for comparison to INLA)
res_idm <- basis.search(base_idm, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10, max.basis.functions = 200)
res_idm <- idm_scampr_w_prediction(dat_pa = dat_pa, dat_po = dat_po, domain = dat)

# figure to compare coefficient estimates requires standard errors
intervals_pa <- confint(res_pa)
intervals_po <- confint(res_po)
intervals_idm <- confint(res_idm)
# get the common term names
fixed.terms <- intersect(names(res_pa$coefficients), names(res_idm$coefficients))
# remove the prior variance
fixed.terms <- fixed.terms[fixed.terms != "Prior log sd(u) (res. 1)"]
intervals_pa <- intervals_pa[fixed.terms, ]
intervals_po <- intervals_po[fixed.terms, ]
intervals_idm <- intervals_idm[fixed.terms, ]

intA <- data.frame(setNames( data.frame(intervals_pa),
                             c("lo", "hi")), Model = "A", est = res_pa$coefficients[names(res_pa$coefficients) %in% fixed.terms], var = fixed.terms)
intB <- data.frame(setNames( data.frame(intervals_po),
                             c("lo", "hi")), Model = "B", est = res_po$coefficients[names(res_po$coefficients) %in% fixed.terms], var = fixed.terms)
intC <- data.frame(setNames( data.frame(intervals_idm),
                             c("lo", "hi")), Model = "C", est = res_idm$coefficients[names(res_idm$coefficients) %in% fixed.terms], var = fixed.terms)
plot_dat <- rbind(intA, intB, intC)


# or
intervals_pa <- data.frame(res_pa$fixed.effects[fixed.terms, ])
intervals_po <- data.frame(res_po$fixed.effects[fixed.terms, ])
intervals_idm <- data.frame(res_idm$fixed.effects[fixed.terms, ])
colnames(intervals_pa) <- c("est", "se")
colnames(intervals_po) <- c("est", "se")
colnames(intervals_idm) <- c("est", "se")
intervals_pa$lo <- intervals_pa$est - intervals_pa$se
intervals_po$lo <- intervals_po$est - intervals_po$se
intervals_idm$lo <- intervals_idm$est - intervals_idm$se
intervals_pa$hi <- intervals_pa$est + intervals_pa$se
intervals_po$hi <- intervals_po$est + intervals_po$se
intervals_idm$hi <- intervals_idm$est + intervals_idm$se

intA <- data.frame(intervals_pa, Model = "A", est = res_pa$coefficients[names(res_pa$coefficients) %in% fixed.terms], var = fixed.terms)
intB <- data.frame(intervals_po, Model = "B", est = res_po$coefficients[names(res_po$coefficients) %in% fixed.terms], var = fixed.terms)
intC <- data.frame(intervals_idm, Model = "C", est = res_idm$coefficients[names(res_idm$coefficients) %in% fixed.terms], var = fixed.terms)
plot_dat <- rbind(intA, intB, intC)

# factor the names so they are in a decent order
plot_dat$var <- factor(plot_dat$var, levels = fixed.terms)

# Plotting the intervals
plot.res <- 500
png(filename = paste0(home.wd, "/app_estimates.png"), width = 6*plot.res, height = 5*plot.res, res = plot.res)
ggplot(plot_dat %>% filter(plot_dat$var != "(Intercept)" & plot_dat$Model != "B" & !plot_dat$var %in% c("tempann", "I(tempann^2)")), aes(x = var, y = est, color = Model)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  labs(title = paste0("Coefficient Estimates for ", sp.numbers$full_name[sp.numbers$job == job]),
       x = "Predictors",
       y = "Estimate (\u00B1 Standard Errors)",
       color = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_blank()) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_color_manual(values = c("coral", "royalblue1"), labels = c("PA only", "IDM"))
dev.off()

# INLA
bnd.time <- system.time(assign("nsw_bnd", inla.nonconvex.hull(as.matrix(dat[ , c("x", "y")]), convex = 0.05)))
mesh <- inla.mesh.2d(loc.domain = dat[ , c("x", "y")], max.edge=c(10,30), cutoff=2, offset = c(5,20))
mesh <- inla.mesh.2d(loc.domain = dat[ , c("x", "y")], max.edge=c(0.5,1), cutoff=2, offset = c(0.2,0.5))
mesh <- inla.mesh.2d(loc.domain = dat[ , c("x", "y")], max.edge=c(0.25,0.5), cutoff=2, offset = c(0.5,1))
mesh <- inla.mesh.2d(loc.domain = dat[ , c("x", "y")], max.edge=c(0.25,0.5), cutoff=0.01, offset = c(0.5,1)) # n.int = 324
mesh <- inla.mesh.2d(loc.domain = dat[ , c("x", "y")], max.edge=c(0.25,0.5), cutoff=0.05)
mesh.time <- system.time(assign("mesh", inla.mesh.2d(loc.domain = dat[ , c("x", "y")], max.edge=c(0.35,0.5), cutoff=0.05)))

check.interior <- mgcv::in.out(as.matrix(nsw_bnd$loc), as.matrix(mesh$loc[,1:2]))
sum(check.interior)
plot(mesh)
lines(nsw_bnd, col = "red", lwd = 2)
points(mesh$loc[,1:2], col = check.interior + 1, pch = c(1, 16)[check.interior + 1])

# set the spde representation to be the mesh with defaults for parameters
spde <- inla.spde2.matern(mesh)
# make A matrix for structured data
structured_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(dat_pa[ , c("x","y")]))

# make A matrix for unstructured data
unstructured_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(dat_po[dat_po$occ == 1 , c("x","y")]))

# make A matrix for quadrature points
quad_data_A <- inla.spde.make.A(mesh = mesh, loc = as.matrix(dat_po[dat_po$occ == 0 , c("x","y")]))

# set the number of integration points, presence points in PO data and prediction points
nq <- sum(dat_po$occ == 0)
n <- sum(dat_po$occ == 1)

# change data to include 0s for nodes and 1s for presences
y.pp <- dat_po$occ

# add expectation vector (area for integration points/nodes and 0 for presences)
e.pp <- dat_po$quad.size

# make the full A matrix for the PPM component
A.pp <- rbind(quad_data_A, unstructured_data_A)

# need to make our own factor contrasts for INLA
fixed.po <- data.frame(get.design.matrix(res_idm$formula, dat_po))
fixed.pa <- data.frame(get.design.matrix(res_idm$formula, dat_pa))

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

# combine the stacks
stk <- inla.stack(stk_unstructured_data, stk_structured_data)

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


ggplot(plot_dat %>% filter(plot_dat$var != "(Intercept)" & plot_dat$Model != "B" & !plot_dat$var %in% c("tempann", "I(tempann^2)")), aes(x = Model, y = est, color = Model)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) + facet_wrap(~ var, scales = "free_y")
  labs(title = "Comparison of Coefficient Estimates",
       x = "Variables",
       y = "Confidence Intervals",
       color = "Models") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_blank()) +
  geom_hline(yintercept = 0, lty = "dashed")

