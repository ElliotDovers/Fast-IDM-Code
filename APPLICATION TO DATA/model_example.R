home.wd <- getwd()
dir.create(file.path(home.wd, "figures"))

plot.res <- 500

# get functions to fit IDMs using scampr and INLA and make predictions across the domain
source("idm_inla.R")
source("idm_scampr.R")
source("idm_mgcv.R")

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
################################################################################

# Perform the spatial k-fold cross-validation #

# Here we are interested in demonstrating the benefit of the IDM for a single species
# job = 7 # Eucalyptus campanulata
# job = 2 # Corymbia gummifera
job = 6 # Eucalyptus fastigata

# get the data in the appropriate format
source("get_and_format_flora_data.R")

## Model the particular species ##

# obtain the species
s = species[job]

# subset presence records of species for this species
presences <- presences[presences$spid == s, ]

# add in the quadrature weights for the background pts (presence weights will be set in fitting routines)
background$quad.size <- region_size / nrow(background)
domain$quad.size <- region_size / nrow(domain)
dat_pa$quad.size <- 1

# combine for the mgcv/scampr friendly data frame
presences$quad.size <- 1e-6
dat_po <- rbind(presences, background)

# subsets the PA dataset to this particular species
dat_pa <- dat_pa[dat_pa$spid == s, ]

# retain the original predictor set
og_preds <- preds

### Fit the models #############################################################

# remove highly correlated variables
library(corrplot)
nsw_cor <- cor(domain[ , preds[!preds %in% c("disturb", "soilfert")]])
corrplot(nsw_cor)

# removing Moisture Index and rainfall in driest quarter. Also removing tempmin as this is highly correlated with tempann
rm.preds <- c("mi", "raindq", "tempmin")
preds <- og_preds[!og_preds %in% rm.preds]

# perform forward selection on the fixed effects
pa.glm0 <- glm(occ ~ 1,
               family = binomial(link = "cloglog"), data = dat_pa
)
pa.forward <- stepAIC(pa.glm0, scope = as.formula(paste0("occ ~ ", paste(preds, collapse = " + "), " + ", paste(paste0("I(", preds[!preds %in% c("disturb", "soilfert")], "^2)"), collapse = " + "))),
                      trace = F, direction = "forward"
)
form <- pa.forward$formula
glm_pa <- pa.forward
# form <- occ ~ tempann + rainann + I(tempann^2) + I(rainann^2) + soilfert
# glm_pa <- glm(form, family = binomial(link = "cloglog"), data = dat_pa)

library(statmod)
set.seed(1)
res_pa = qresiduals(glm_pa)
set.seed(NULL)

# add in a GRF for spatial error
library(mgcv)
sp.glm_pa <- gam(update(form, . ~ . + s(x,y,bs="gp",k=200)), family = binomial(link = "cloglog"), data = dat_pa, method = "REML")
set.seed(1)
sp.res_pa = qresiduals(sp.glm_pa)
set.seed(NULL)

# check the residual correlograms
library(ncf)
if (file.exists("E_fast_correlogs.RDATA")) {
  load("E_fast_correlogs.RDATA")
} else {
  cl_pa = spline.correlog(x = dat_pa$x, y = dat_pa$y, z = res_pa)
  sp.cl_pa = spline.correlog(x = dat_pa$x, y = dat_pa$y, z = sp.res_pa)
  save(list = c("cl_pa", "sp.cl_pa"), file = "E_fast_correlogs.RDATA")
}

# plot the check for spatial correlation
png(filename = paste0(home.wd, "/figures/correlogs.png"), width = 6*plot.res, height = 4*plot.res, res = plot.res)
par(mfrow = c(2, 1), mar = c(1.5,3.1,1.6,0))
plot(cl_pa, ylim = c(-0.15,0.15), xaxt = "n", yaxt = "n", xlim = c(0,100))
axis(2, at = c(-0.1, 0, 0.1))
mtext("Correlation", 2, line = 2.3)
mtext("A: model with predictors only", 3, adj = 0)
par(mar = c(3.1,3.1,0,0))
plot(sp.cl_pa, ylim = c(-0.15,0.15), yaxt = "n", xlim = c(0,100))
axis(2, at = c(-0.1, 0, 0.1))
mtext("Correlation", 2, line = 2.3)
mtext("B: model with predictors and a latent spatial field", 3, adj = 0)
mtext("distance (km)", 1, line = 2.2)
par(mfrow = c(1, 1))
dev.off()

# Fit the ISDMs for comparison of software

if (file.exists("applied_example.RDATA")) {
  load("applied_example.RDATA")
} else {
  # fit the IDM via scampr - within a function for timing (note that for scampr this fits base_idm again to include the comp. time)
  scampr.time <- system.time(assign("res_idm_scampr", idm_scampr(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain)))
  # save(res_idm_scampr, file = "res_idm_scampr.RDATA")
  # remove unecessary env objects to save space
  res_idm_scampr$call <- NULL
  res_idm_scampr$tmb.call.list <- NULL
  attr(res_idm_scampr$formula, ".Environment") <- NULL
  attr(attr(res_idm_scampr$formula, "bias"), ".Environment") <- NULL
  # bfs <- simple_basis(sqrt(300), domain)
  # scampr.time <- system.time(assign("res_idm_scampr", scampr(form, data = dat_po, bias.formula = ~ 1, pa.data = dat_pa, basis.functions = bfs, latent.po.biasing = T, model.type = "IDM")))
  
  # fit the IDM via mgcv - within a function for timing
  mgcv.time <- system.time(assign("res_idm_mgcv", idm_mgcv(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain, fld_dim = 200, bias_fld_dim = 200, estimate.range = T)))
  # save(res_idm_mgcv, file = "res_idm_mgcv.RDATA")
  
  # fit the IDM via INLA - within a function for timing
  inla.time <- system.time(assign("res_idm_inla", idm_inla(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain, fld_dim = 200)))
  # save(res_idm_inla, file = "res_idm_inla.RDATA")
  
  save(list = c("res_idm_scampr", "res_idm_mgcv", "res_idm_inla", "scampr.time", "mgcv.time", "inla.time"), file = "applied_example.RDATA")
}

# get the common term names for the fixed effect comparisons
fixed.terms <- intersect(names(res_idm_scampr$coefficients), names(res_idm_mgcv$coefficients))
fixed.terms <- fixed.terms[fixed.terms != "(Intercept)"]
fixed.terms_inla <- gsub("^", ".", gsub(")", ".", gsub("(", ".", fixed.terms, fixed = T), fixed = T), fixed = T)

# get the interval estimates
inla_ints <- cbind(res_idm_inla$summary.fixed[fixed.terms_inla, "mean"] + qnorm(0.025) * res_idm_inla$summary.fixed[fixed.terms_inla, "sd"], res_idm_inla$summary.fixed[fixed.terms_inla, "mean"] + qnorm(0.975) * res_idm_inla$summary.fixed[fixed.terms_inla, "sd"])
intA <- data.frame(lo = inla_ints[ , 1], hi = inla_ints[ , 2], fit = "R-INLA", est = res_idm_inla$summary.fixed[fixed.terms_inla, "mean"], var = fixed.terms)

scampr_ints <- confint(res_idm_scampr)
intB <- data.frame(lo = scampr_ints[fixed.terms, 1], hi = scampr_ints[fixed.terms, 2], fit = "scampr", est = coef(res_idm_scampr)[fixed.terms], var = fixed.terms, row.names = 1:length(fixed.terms))

tmp <- summary(res_idm_mgcv)
mgcv_ints <- cbind(coef(res_idm_mgcv)[fixed.terms] + qnorm(0.025) * tmp$se[fixed.terms], coef(res_idm_mgcv)[fixed.terms] + qnorm(0.975) * tmp$se[fixed.terms])
intC <- data.frame(lo = mgcv_ints[ , 1], hi = mgcv_ints[ , 2], fit = "mgcv", est = coef(res_idm_mgcv)[fixed.terms], var = fixed.terms)
plot_dat <- rbind(intA, intB, intC)

# factor the names so they are in a decent order
plot_dat$var <- factor(plot_dat$var, levels = rev(fixed.terms)[c(2,1,3:length(fixed.terms))])
plot_dat$fit <- factor(plot_dat$fit, levels = c("R-INLA", "scampr", "mgcv"))

# Plotting the intervals
library(ggplot2)
plot.res <- 500
png(filename = paste0(home.wd, "/figures/app_estimates.png"), width = 6*plot.res, height = 5*plot.res, res = plot.res)
ggplot(plot_dat, aes(x = var, y = est, color = fit)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.5, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5)) +
  labs(title = "Comparison of Wald confidence intervals (\u03B1 = 0.05)",
       x = "",
       y = "Estimate",
       color = "Software") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_color_manual(values = c("darkorange1", "dodgerblue3", "darkorchid4")) + coord_flip()
dev.off()

# predict the relative abundance surfaces from each model fit:

# add necessary elements to domain
domain$source <- 1
domain$po_id <- 0

fld_mgcv <- predict(res_idm_mgcv, newdata = domain)
fld_scampr <- predict(res_idm_scampr, newdata = domain, include.bias.accounting = F)
# for inla we need to re-fit the model
tmp <- idm_inla(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain, fld_dim = 200, with.predictions = T)
fld_inla <- log(tmp$preds)

# set some things for the plot
library(viridis)

# truncate the estimates to -10 (there are only a small location in the south-west mountainous region for each model that gets values this low)
fld_mgcv[fld_mgcv < -10] <- -10
fld_scampr[fld_scampr < -10] <- -10
fld_inla[fld_inla < -10] <- -10

zlims <- range(fld_mgcv, fld_scampr, fld_inla)
ribbon.col.width <- 0.15
png(filename = paste0(home.wd, "/figures/fld_estimates.png"), width = 6*plot.res, height = 3.5*plot.res, res = plot.res)
layout(matrix(1:4, ncol = 4), widths = c(rep((1 - ribbon.col.width) / 3, 3), ribbon.col.width))
# par(mfrow = c(1, 3), mar = c(8,0,8,0))
par(mar = c(0,0,2,0))
plot(vec2im(fld_inla, domain$x, domain$y), box = F, main = "", zlim = zlims, ribbon = F, col = inferno)
title(main = "A: R-INLA")
plot(vec2im(fld_scampr, domain$x, domain$y), box = F, main = "", zlim = zlims, ribbon = F, col = inferno)
title(main = "B: scampr")
plot(vec2im(fld_mgcv, domain$x, domain$y), box = F, main = "", zlim = zlims, ribbon = F, col = inferno)
title(main = "C: mgcv")
par(mar = c(0,2,4,1))
legend_image <- as.raster(matrix(rev(inferno(20)), ncol=1))
plot(c(0,2),c(zlims[1],zlims[2]),type = 'n', axes = F, xlab = "", ylab = "", main = "")
mtext(expression(log~hat(mu)(bold(s))), line = 0, cex = 0.75, adj = 0)
text(x=1.75, y = seq(zlims[1],0,l=6), labels = seq(-10,0,l=6))
lines(c(0,1.01,1.01,0,0), c(zlims[1],zlims[1],zlims[2],zlims[2],zlims[1]))
for (i in seq(-10,0,l=6)) {
  lines(c(1.01, 1.21), rep(i, 2))
}
rasterImage(legend_image, 0, zlims[1], 1, zlims[2])
par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()

# also obtain the latent fields
domain$po_id <- 1
tmp <- predict(res_idm_mgcv, newdata = domain, type = "terms")
lat_fld_mgcv <- tmp[,"s(x,y)"]
bias_fld_mgcv <- tmp[,"s(x,y):po_id"]
lat_fld_scampr <- attr(fld_scampr, "Zmu")

library(viridis)

# Plot up some of the fields for visualising the ISDM:

png(filename = paste0(home.wd, "/figures/pa.png"), width = 5*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(domain$source, domain$x, domain$y), box = F, main = "", col = "grey", ribbon = F)
points(dat_pa[dat_pa$occ == 0, c("x", "y")], pch = 1, col = "goldenrod")
points(dat_pa[dat_pa$occ == 1, c("x", "y")], pch = 4, col = "royalblue")
legend(x = 180, y = 6920, legend = c("present", "absent"), pch = c(4,1), col = c("royalblue", "goldenrod"), bty = "n", horiz = F, cex = 3, xpd = T, pt.lwd = 2, x.intersp = c(0.5, 0.5), text.width = c(.23,.23))
dev.off()

png(filename = paste0(home.wd, "/figures/po.png"), width = 5*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(domain$source, domain$x, domain$y), box = F, main = "", ribbon = F, col = "grey")
points(presences[presences$spid == species[job], c("x", "y")], pch = 4, col = "magenta4")
legend(x = 175, y = 6920, legend = c("presence\nonly"), pch = 4, col = "magenta4", bty = "n", horiz = F, cex = 3, xpd = T, pt.lwd = 2, x.intersp = c(0.5), text.width = c(.23))
dev.off()

png(filename = paste0(home.wd, "/figures/mu.png"), width = 4.2*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(fld_mgcv, domain$x, domain$y), box = F, main = "", ribbon = F, col = inferno)
# text(x = 320, y = 6820, label = expression(paste("ln", mu(bold(s)))), cex = 3.5, xpd = T)
text(x = 320, y = 6820, label = "Mean", cex = 3.5, xpd = T)
dev.off()

png(filename = paste0(home.wd, "/figures/xi.png"), width = 4.2*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(lat_fld_mgcv, domain$x, domain$y), box = F, main = "", ribbon = F, col = plasma)
# text(x = 320, y = 6820, label = expression(xi(bold(s))), cex = 3.5, xpd = T)
text(x = 320, y = 6820, label = "Error", cex = 3.5, xpd = T)
dev.off()

png(filename = paste0(home.wd, "/figures/xi_B.png"), width = 4.2*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(bias_fld_mgcv, domain$x, domain$y), box = F, main = "", ribbon = F, col = cividis)
# text(x = 320, y = 6820, label = expression(xi[B](bold(s))), cex = 3.5, xpd = T)
text(x = 320, y = 6820, label = "Bias", cex = 3.5, xpd = T)
dev.off()

png(filename = paste0(home.wd, "/figures/env.png"), width = 4.2*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(apply(domain[,all.vars(form)[2:length(all.vars(form))]], 1, function(x){sum(as.numeric(x))}), domain$x, domain$y), box = F, main = "", ribbon = F, col = terrain.colors)
# text(x = 320, y = 6820, label = expression(X(bold(s))), cex = 3.5, xpd = T)
text(x = 300, y = 6820, label = "Env.", cex = 3.5, xpd = T)
dev.off()

## Collate the coefficient estimates

mgcv_ests <- summary(res_idm_mgcv)$p.table[-nrow(summary(res_idm_mgcv)$p.table), c("Estimate", "Std. Error")]
scampr_ests <- res_idm_scampr$fixed.effects
inla_ests <- res_idm_inla$summary.fixed[-2,c("mean", "sd")]
row.names(inla_ests)[row.names(inla_ests) == "interceptA"] <- "(Intercept)"

# check the orders
identical(row.names(mgcv_ests), row.names(scampr_ests))
identical(row.names(mgcv_ests), row.names(inla_ests))

tab <- cbind(paste0(apply(round(mgcv_ests, 3), 1, paste, collapse = " ("), ")"),
             paste0(apply(round(scampr_ests, 3), 1, paste, collapse = " ("), ")"),
             paste0(apply(round(inla_ests, 3), 1, paste, collapse = " ("), ")")
)
dimnames(tab) <- list(row.names(mgcv_ests), c("mgcv", "scampr", "R-INLA"))
library(xtable)
print(xtable(tab, caption = "Comparison of fixed effect estimates (with uncertainty represented by standard errors in brackets) from IDMs fitted via \texttt{mgcv}, \texttt{scampr} and \texttt{INLA} modelling occurence data for \textit{Eucalyptus campanulata} in Northern NSW, Australia.", label = "tab:coeff_comp", align = c("l", "r", "r", "r")))




# fit the base scampr models
base_po <- scampr(form, data = dat_po, include.sre = F, model.type = "PO", sre.approx = "laplace")
base_pa <- scampr(form, data = dat_pa, include.sre = F, model.type = "PA", sre.approx = "laplace")
base_idm <- scampr(form, data = dat_po, bias.formula = ~ 1, pa.data = dat_pa, include.sre = F, model.type = "IDM", sre.approx = "laplace", latent.po.biasing = F)

# fit the PA model
m_pa <- basis.search.pa(base_pa, domain.data = domain, return.model = T, start.nodes = 10)
# fit the PO model
m_po <- basis.search.po(base_po, domain.data = domain, return.model = T, start.nodes = 10)

# source("pa_mgcv.R")
# source("po_mgcv.R")
# 
# # fit the PA Only model via mgcv - within a function for timing
# inla_pa.time <- system.time(assign("res_pa_mgcv", pa_mgcv(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain, fld_dim = 100)))
# 
# # fit the PO Only model via mgcv - within a function for timing
# mgcv_po.time <- system.time(assign("res_po_mgcv", po_mgcv(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain, fld_dim = 100)))

# set the model formula
form <- as.formula(paste0("occ ~ ", paste(preds, collapse = " + "), " + ", paste(paste0("I(", preds[!preds %in% c("disturb", "soilfert")], "^2)"), collapse = " + ")))

if (file.exists("applied_example.RDATA")) {
  load("applied_example.RDATA")
} else {
  # fit the IDM via scampr - within a function for timing (note that for scampr this fits base_idm again to include the comp. time)
  scampr.time <- system.time(assign("res_idm_scampr", idm_scampr(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain)))
  # save(res_idm_scampr, file = "res_idm_scampr.RDATA")
  
  # fit the IDM via mgcv - within a function for timing
  mgcv.time <- system.time(assign("res_idm_mgcv", idm_mgcv(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain, fld_dim = 200)))
  mgcv.time <- system.time(assign("res_idm_mgcv", idm_mgcv(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain, fld_dim = 112, bias_fld_dim = 27)))
  # save(res_idm_mgcv, file = "res_idm_mgcv.RDATA")
  
  # fit the IDM via INLA - within a function for timing
  inla.time <- system.time(assign("res_idm_inla", idm_inla(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain, fld_dim = 200)))
  # save(res_idm_inla, file = "res_idm_inla.RDATA")
  
  save(list = c("res_idm_scampr", "res_idm_mgcv", "res_idm_inla", "scampr.time", "mgcv.time", "inla.time"), file = "applied_example.RDATA")
}

## Collate the coefficient estimates

mgcv_ests <- summary(res_idm_mgcv)$p.table[-nrow(summary(res_idm_mgcv)$p.table), c("Estimate", "Std. Error")]
scampr_ests <- res_idm_scampr$fixed.effects
inla_ests <- res_idm_inla$summary.fixed[-2,c("mean", "sd")]
row.names(inla_ests)[row.names(inla_ests) == "interceptA"] <- "(Intercept)"

# check the orders
identical(row.names(mgcv_ests), row.names(scampr_ests))
identical(row.names(mgcv_ests), row.names(inla_ests))

tab <- cbind(paste0(apply(round(mgcv_ests, 3), 1, paste, collapse = " ("), ")"),
             paste0(apply(round(scampr_ests, 3), 1, paste, collapse = " ("), ")"),
             paste0(apply(round(inla_ests, 3), 1, paste, collapse = " ("), ")")
)
dimnames(tab) <- list(row.names(mgcv_ests), c("mgcv", "scampr", "R-INLA"))
library(xtable)
print(xtable(tab, caption = "Comparison of fixed effect estimates (with uncertainty represented by standard errors in brackets) from IDMs fitted via \texttt{mgcv}, \texttt{scampr} and \texttt{INLA} modelling occurence data for \textit{Eucalyptus campanulata} in Northern NSW, Australia.", label = "tab:coeff_comp", align = c("l", "r", "r", "r")))

# plot of the parameter estimates
fill_cols <- c("darkorange1", "dodgerblue3", "aquamarine4", "darkorchid4")
plot.res <- 500

png(filename = paste0(home.wd, "/comp_estimates.png"), width = 6*plot.res, height = 3.5*plot.res, res = plot.res)
par(mfrow = c(1, 2), mar = c(3,4,0.5,0))
plot(scampr_ests[,1], mgcv_ests[,1], ylab = "mgcv", xlab = "", asp = 1, pch = ".",
     ylim = range(c(mgcv_ests[,1] - mgcv_ests[,2], mgcv_ests[,1] + mgcv_ests[,2])),
     xlim = range(c(scampr_ests[,1] - scampr_ests[,2], scampr_ests[,1] + scampr_ests[,2])))
title(xlab = "scampr", line = 1.9)
arrows(y0 = mgcv_ests[,1] - mgcv_ests[,2], y1 = mgcv_ests[,1] + mgcv_ests[,2], x0 = scampr_ests[,1], x1 = scampr_ests[,1],
       code = 3, angle = 90, length = 0.05, col = "dodgerblue3")
arrows(x0 = scampr_ests[,1] - scampr_ests[,2], x1 = scampr_ests[,1] + scampr_ests[,2], y0 = mgcv_ests[,1], y1 = mgcv_ests[,1],
       code = 3, angle = 90, length = 0.05, col = "darkorange1")
abline(0,1, lty = "dashed")
par(mar = c(3,2,0.5,2))
plot(inla_ests[,1], mgcv_ests[,1], xlab = "", yaxt = "n", ylab = "", asp = 1, pch = ".",
     ylim = range(c(mgcv_ests[,1] - mgcv_ests[,2], mgcv_ests[,1] + mgcv_ests[,2])),
     xlim = range(c(inla_ests[,1] - inla_ests[,2], inla_ests[,1] + inla_ests[,2])))
title(xlab = "R-INLA", line = 1.9)
arrows(y0 = mgcv_ests[,1] - mgcv_ests[,2], y1 = mgcv_ests[,1] + mgcv_ests[,2], x0 = inla_ests[,1], x1 = inla_ests[,1],
       code = 3, angle = 90, length = 0.05, col = "dodgerblue3")
arrows(x0 = inla_ests[,1] - inla_ests[,2], x1 = inla_ests[,1] + inla_ests[,2], y0 = mgcv_ests[,1], y1 = mgcv_ests[,1],
       code = 3, angle = 90, length = 0.05, col = "darkorange1")
abline(0,1, lty = "dashed")
par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()

# predict the relative abundance surfaces from each model fit:

# add necessary elements to domain
domain$source <- 1
domain$po_id <- 0

fld_mgcv <- predict(res_idm_mgcv, newdata = domain)
fld_scampr <- predict(res_idm_scampr, newdata = domain, include.bias.accounting = F)
# for inla we need to re-fit the model
tmp <- idm_inla(form = form, dat_pa = dat_pa, pres = presences, quad = background, domain = domain, fld_dim = 200, with.predictions = T)
fld_inla <- log(tmp$preds)

# truncate the predictions as there are large negative values on the log-scale throwing things off visually
fld_mgcv[fld_mgcv < -10] <- -10
fld_scampr[fld_scampr < -10] <- -10
fld_inla[fld_inla < -10] <- -10

zlims <- range(fld_mgcv, fld_scampr, fld_inla)

# get the latent fields
domain$po_id <- 1
tmp <- predict(res_idm_mgcv, newdata = domain, type = "terms")
lat_fld_mgcv <- tmp[,"s(x,y)"]
bias_fld_mgcv <- tmp[,"s(x,y):po_id"]
lat_fld_scampr <- attr(fld_scampr, "Zmu")

library(viridis)

plot.res <- 500
ribbon.col.width <- 0.15
png(filename = paste0(home.wd, "/fld_estimates.png"), width = 6*plot.res, height = 3.5*plot.res, res = plot.res)
layout(matrix(1:4, ncol = 4), widths = c(rep((1 - ribbon.col.width) / 3, 3), ribbon.col.width))
# par(mfrow = c(1, 3), mar = c(8,0,8,0))
par(mar = c(0,0,2,0))
plot(vec2im(fld_mgcv, domain$x, domain$y), box = F, main = "", zlim = zlims, ribbon = F, col = inferno)#, ribside = "bottom")
# mtext("A: mgcv", side = 3, line = 1.5, adj = 1, padj = 1)
title(main = "A: mgcv")
plot(vec2im(fld_scampr, domain$x, domain$y), box = F, main = "", zlim = zlims, ribbon = F, col = inferno)#, ribside = "bottom")
# mtext("B: scampr", side = 3, line = 1.5, adj = 1, padj = 1)
title(main = "B: scampr")
plot(vec2im(fld_inla, domain$x, domain$y), box = F, main = "", zlim = zlims, ribbon = F, col = inferno)#, ribside = "bottom")
# mtext("C: R-INLA", side = 3, line = 1.5, adj = 1, padj = 1)
title(main = "C: R-INLA")
par(mar = c(0,2,4,1))
legend_image <- as.raster(matrix(rev(inferno(20)), ncol=1))
plot(c(0,2),c(zlims[1],zlims[2]),type = 'n', axes = F, xlab = "", ylab = "", main = "")
mtext(expression(log~hat(mu)(bold(s))), line = 0, cex = 0.75, adj = 0)
text(x=1.75, y = seq(zlims[1],0,l=6), labels = seq(-10,0,l=6))
lines(c(0,1.01,1.01,0,0), c(zlims[1],zlims[1],zlims[2],zlims[2],zlims[1]))
for (i in seq(-10,0,l=6)) {
  lines(c(1.01, 1.21), rep(i, 2))
}
rasterImage(legend_image, 0, zlims[1], 1, zlims[2])
par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()

## For ESA2024

library(oz)
library(viridis)
library(scampr)

# tree <- dat_pa[dat_pa$spid == species[job], ]
# save(list = c("domain", "tree", "presences", "species.info", "fld_mgcv", "fld_scampr", "fld_inla", "lat_fld_mgcv", "bias_fld_mgcv", "zlims", "background", "species", "job"), file = "esa.RDATA")
load("esa.RDATA")

home.wd <- getwd()
plot.res <- 500
ribbon.col.width <- 0.15

png(filename = paste0(home.wd, "/esa_mean_fld.png"), width = 4*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(fld_mgcv, domain$x, domain$y), box = F, main = "", ribbon = F, col = inferno)
dev.off()

png(filename = paste0(home.wd, "/esa_pa.png"), width = 4.2*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(domain$source, domain$x, domain$y), box = F, main = "", ribbon = F, col = "grey")
points(tree[tree$occ == 0, c("x", "y")], pch = 1, col = "goldenrod")
points(tree[tree$occ == 1, c("x", "y")], pch = 4, col = "royalblue")
legend(x = 220, y = 6865, legend = c("presence", "absence"), pch = c(4,1), col = c("royalblue", "goldenrod"), bty = "n", horiz = F, cex = 2, xpd = T)
dev.off()

png(filename = paste0(home.wd, "/esa_po.png"), width = 4.2*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(domain$source, domain$x, domain$y), box = F, main = "", ribbon = F, col = "grey")
points(presences[presences$spid == species[job], c("x", "y")], pch = 4, col = "magenta4")
legend(x = 220, y = 6865, legend = c("presence\nonly"), pch = 4, col = "magenta4", bty = "n", horiz = F, cex = 2, xpd = T)
dev.off()

png(filename = paste0(home.wd, "/esa_nsw.png"), width = 6*plot.res, height = 4*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
oz()
plot(scampr::vec2im(domain$source, domain$lon, domain$lat), box = F, ribbon = F, main = "", add = T, col = "grey")
points(y = -30.2749, x = 153.1337, pch = NA, cex = 3, col = "red", lwd = 5)
dev.off()

png(filename = paste0(home.wd, "/esa_bias_fld.png"), width = 4*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(bias_fld_mgcv, domain$x, domain$y), box = F, main = "", ribbon = F, col = cividis)
dev.off()

png(filename = paste0(home.wd, "/esa_climate.png"), width = 4*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(-(fld_mgcv - lat_fld_mgcv), domain$x, domain$y), box = F, main = "", ribbon = F, col = terrain.colors)
# points(presences[presences$spid == species[job], c("x", "y")], pch = 4, col = "magenta4")
# legend(x = 220, y = 6865, legend = c("presence\nonly"), pch = 4, col = "magenta4", bty = "n", horiz = F, cex = 2, xpd = T)
dev.off()

png(filename = paste0(home.wd, "/esa_error.png"), width = 4*plot.res, height = 6*plot.res, res = plot.res)
par(mar = c(0,0,0,0))
plot(vec2im(lat_fld_mgcv, domain$x, domain$y), box = F, main = "", ribbon = F, col = plasma)
dev.off()

png(filename = paste0(home.wd, "/esa_mean_comp.png"), width = 6*plot.res, height = 3.5*plot.res, res = plot.res)
layout(matrix(1:4, ncol = 4), widths = c(rep((1 - ribbon.col.width) / 3, 3), ribbon.col.width))

par(mar = c(0,0,2,0))
plot(vec2im(fld_inla, domain$x, domain$y), box = F, main = "", zlim = zlims, ribbon = F, col = inferno)
# title(main = "R-INLA")

plot(vec2im(fld_scampr, domain$x, domain$y), box = F, main = "", zlim = zlims, ribbon = F, col = inferno)
# title(main = "scampr")

plot(vec2im(fld_mgcv, domain$x, domain$y), box = F, main = "", zlim = zlims, ribbon = F, col = inferno)
# title(main = "mgcv")


par(mar = c(0,2,4,1))
legend_image <- as.raster(matrix(rev(inferno(20)), ncol=1))
plot(c(0,2),c(zlims[1],zlims[2]),type = 'n', axes = F, xlab = "", ylab = "", main = "")
# mtext(expression(log~"Underlying Mean"), line = 0, cex = 0.75, adj = 0)
text(x=1.75, y = seq(zlims[1],0,l=6), labels = seq(-10,0,l=6))
lines(c(0,1.01,1.01,0,0), c(zlims[1],zlims[1],zlims[2],zlims[2],zlims[1]))
for (i in seq(-10,0,l=6)) {
  lines(c(1.01, 1.21), rep(i, 2))
}
rasterImage(legend_image, 0, zlims[1], 1, zlims[2])
par(mfrow = c(1, 1), mar = c(5.1,4.1,4.1,2.1))
dev.off()


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