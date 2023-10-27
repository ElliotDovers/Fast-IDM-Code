home.wd <- getwd()

# get functions to fit IDMs using scampr and INLA and make predictions across the domain
source("idm_inla.R")
source("idm_scampr.R")

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
job = 7 # Eucalyptus campanulata

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

# load up the pre-calculated means and standard deviations of predictors (across the entire domain) for scaling
load("scaling_info.RDATA")
load("nsw_grid.RDATA") # grid points of the entire domain

# # get the full domain data grid
# load("NSW.RDATA") # load the pre-prepared full domain data
# covars <- preds[!preds %in% c("disturb", "soilfert")]
# covar_means <- colMeans(dat[ , covars], na.rm = T)
# covar_sds <- sapply(dat[ , covars], sd, na.rm = T)
# # calculate the scaled variables
# tmp_covars_scaled <- scale(dat[, covars], center = covar_means, scale = covar_sds)
# # retain the original variables
# colnames(dat)[colnames(dat) %in% covars] <- paste0(colnames(dat)[colnames(dat) %in% covars], "_og")
# # add in the scaled variables
# dat <- cbind(dat, tmp_covars_scaled)
# rm(tmp_covars_scaled)
# gc()

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
      # dat[ ,fac_col][dat[ ,fac_col] %in% c(4,5)] <- 3
    }
    # dat[ ,fac_col] <- as.factor(dat[ ,fac_col])
    dat_po[ ,fac_col] <- as.factor(dat_po[ ,fac_col])
    dat_pa[ ,fac_col] <- as.factor(dat_pa[ ,fac_col])
    # expand and relevel incase the datasets are missing levels from the full domain
    # dat_po[ ,fac_col] <- forcats::fct_expand(dat_po[,fac_col], levels(dat[,fac_col]))
    # dat_po[ ,fac_col] <- forcats::fct_relevel(dat_po[,fac_col], levels(dat[,fac_col]))
    dat_pa[ ,fac_col] <- forcats::fct_expand(dat_pa[,fac_col], levels(dat_po[,fac_col]))
    dat_pa[ ,fac_col] <- forcats::fct_relevel(dat_pa[,fac_col], levels(dat_po[,fac_col]))
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
base_idm <- scampr(form, data = dat_po, bias.formula = ~ 1, pa.data = dat_pa, include.sre = F, model.type = "IDM", sre.approx = "laplace", latent.po.biasing = F)

# fit the PA model
res_pa <- basis.search.pa(base_pa, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10, max.basis.functions = 200)
# fit the PO model
res_po <- basis.search.po(base_po, domain.data = dat_po[dat_po$occ == 0, ], return.model = T, start.nodes = 10, which.approx = "laplace", max.basis.functions = 100)

# fit the IDM - within a function for timing (note that for scampr this fits base_idm again to include the comp. time)
scampr.time <- system.time(assign("res_idm", idm_scampr(form = form, dat_pa = dat_pa, dat_po = dat_po, domain = dat)))

# also fit the INLA equivalent INLA model
inla.time <- system.time(assign("res_idm_inla", idm_inla(form = form, dat_pa = dat_pa, dat_po = dat_po, domain = dat)))

save(list = c("res_pa", "res_po", "res_idm", "res_idm_inla", "scampr.time", "inla.time"), file = "applied_example.RDATA")

tab <- cbind(paste0(apply(round(res_pa$fixed.effects, 2), 1, paste, collapse = " ("), ")"),
      paste0(apply(round(res_idm$fixed.effects, 2), 1, paste, collapse = " ("), ")"),
      paste0(apply(round(res_idm_inla$summary.fixed[-2,c("mean", "sd")], 2), 1, paste, collapse = " ("), ")")
)
dimnames(tab) <- list(rownames(res_pa$fixed.effects), c("PA only", "IDM - scampr", "IDM - INLA"))
library(xtable)
print(xtable(tab, caption = "Comparison of fixed effect estimates (with uncertainty represented by standard errors in brackets) from a PA only model and IDM  (fitted via \texttt{scampr} and \texttt{INLA}) modelling \textit{Eucalyptus campanulata} in Northern NSW, Australia.", label = "tab:coeff_comp", align = c("l", "r", "r", "r")))

# # figure to compare coefficient estimates requires standard errors
# intervals_pa <- confint(res_pa)
# intervals_po <- confint(res_po)
# intervals_idm <- confint(res_idm)
# # get the common term names
# fixed.terms <- intersect(names(res_pa$coefficients), names(res_idm$coefficients))
# # remove the prior variance
# fixed.terms <- fixed.terms[fixed.terms != "Prior log sd(u) (res. 1)"]
# intervals_pa <- intervals_pa[fixed.terms, ]
# intervals_po <- intervals_po[fixed.terms, ]
# intervals_idm <- intervals_idm[fixed.terms, ]
# 
# intA <- data.frame(setNames( data.frame(intervals_pa),
#                              c("lo", "hi")), Model = "A", est = res_pa$coefficients[names(res_pa$coefficients) %in% fixed.terms], var = fixed.terms)
# intB <- data.frame(setNames( data.frame(intervals_po),
#                              c("lo", "hi")), Model = "B", est = res_po$coefficients[names(res_po$coefficients) %in% fixed.terms], var = fixed.terms)
# intC <- data.frame(setNames( data.frame(intervals_idm),
#                              c("lo", "hi")), Model = "C", est = res_idm$coefficients[names(res_idm$coefficients) %in% fixed.terms], var = fixed.terms)
# plot_dat <- rbind(intA, intB, intC)
# 
# 
# # or
# intervals_pa <- data.frame(res_pa$fixed.effects[fixed.terms, ])
# intervals_po <- data.frame(res_po$fixed.effects[fixed.terms, ])
# intervals_idm <- data.frame(res_idm$fixed.effects[fixed.terms, ])
# colnames(intervals_pa) <- c("est", "se")
# colnames(intervals_po) <- c("est", "se")
# colnames(intervals_idm) <- c("est", "se")
# intervals_pa$lo <- intervals_pa$est - intervals_pa$se
# intervals_po$lo <- intervals_po$est - intervals_po$se
# intervals_idm$lo <- intervals_idm$est - intervals_idm$se
# intervals_pa$hi <- intervals_pa$est + intervals_pa$se
# intervals_po$hi <- intervals_po$est + intervals_po$se
# intervals_idm$hi <- intervals_idm$est + intervals_idm$se
# 
# intA <- data.frame(intervals_pa, Model = "A", est = res_pa$coefficients[names(res_pa$coefficients) %in% fixed.terms], var = fixed.terms)
# intB <- data.frame(intervals_po, Model = "B", est = res_po$coefficients[names(res_po$coefficients) %in% fixed.terms], var = fixed.terms)
# intC <- data.frame(intervals_idm, Model = "C", est = res_idm$coefficients[names(res_idm$coefficients) %in% fixed.terms], var = fixed.terms)
# plot_dat <- rbind(intA, intB, intC)
# 
# # factor the names so they are in a decent order
# plot_dat$var <- factor(plot_dat$var, levels = fixed.terms)
# 
# # Plotting the intervals
# plot.res <- 500
# png(filename = paste0(home.wd, "/app_estimates.png"), width = 6*plot.res, height = 5*plot.res, res = plot.res)
# ggplot(plot_dat %>% filter(plot_dat$var != "(Intercept)" & plot_dat$Model != "B" & !plot_dat$var %in% c("tempann", "I(tempann^2)")), aes(x = var, y = est, color = Model)) +
#   geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.5, position = position_dodge(width = 0.5)) +
#   geom_point(position = position_dodge(width = 0.5)) +
#   labs(title = paste0("Coefficient Estimates for ", sp.numbers$full_name[sp.numbers$job == job]),
#        x = "Predictors",
#        y = "Estimate (\u00B1 Standard Errors)",
#        color = "Model") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_blank()) +
#   geom_hline(yintercept = 0, lty = "dashed") +
#   scale_color_manual(values = c("coral", "royalblue1"), labels = c("PA only", "IDM"))
# dev.off()