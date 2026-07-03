## Collate and Visual Results from the simulation ##

library(dplyr)
library(ggplot2)

# sets which results to compile - the following corresponds to the scenario presented in the main paper
res_type <- "_all"

home.wd <- getwd()
# Get the job array
# tab <- read.csv(paste0("job_array", res_type, ".csv"))

# initialise the result storage
dat <- NULL

res.list <- list()
res.objs <- list.files(paste0(home.wd, "/Results", res_type))[grepl("res_", list.files(paste0(home.wd, "/Results", res_type)), fixed = T)]
for (job in 1:length(res.objs)) {
  load(paste0(home.wd, "/Results", res_type, "/", res.objs[job]))
  res.list[[job]] <- res_tab
}
dat <- do.call(rbind, res.list)
rm(res.list)

# create basis function identifier
dat$BASIS_TYPE <- "Optimised"
dat$BASIS_TYPE[dat$fit_model %in% c("MGCV FIXED", "SCAMPR FIXED")] <- "Fixed"

# adjust the fit to reflect only the software
dat$FIT[dat$FIT == "SCAMPR FIXED"] <- "SCAMPR"
# dat$Model <- factor(paste0(dat$FIT, " ", dat$MODEL),
#                     levels = c("INLA PA", "SCAMPR PA", "SCAMPR2 PA", "INLA IDM", "SCAMPR IDM", "SCAMPR2 IDM", "INLA PO", "SCAMPR PO", "SCAMPR2 PO"))
dat$fit_model[dat$fit_model %in% c("INLA PA", "INLA PO", "INLA IDM")] <- "INLA"
dat$Fit <- factor(dat$fit_model,
                  levels = c("INLA", "SCAMPR FIXED", "SCAMPR", "MGCV FIXED", "MGCV"),
                  labels = c("INLA", "scampr", "scampr (opt. k)","mgcv (def.)", "mgcv (opt.)")
)

# get the rates of failure to converge or other problems
fail.rates.timing <- dat %>% group_by(MODEL, FIT, BASIS_TYPE) %>%
  summarise(poor.conv.kl = sum(!(KL < 1e5 & !is.na(KL))), poor.conv.mae = sum(!(MAE < 1e5 & !is.na(MAE))), cpu = mean(TIME), fails = 100 - length(unique(sim)))
# ploterr <- fail.rates.timing[fail.rates.timing$Model %in% c("INLA IDM", "SCAMPR IDM", "SCAMPR2 IDM"), ]
# ploterr$Model <- factor(ploterr$Model)
# ploterr$fail.rate <- ploterr$fails / 100

ggplot(dat, aes(y = KL, x = Fit)) +
  facet_wrap(~ fld_dim + MODEL, scales = "free") +
  geom_boxplot()

ggplot(dat %>% filter(MODEL == "IDM"), aes(y = KL, x = Fit, fill = Fit)) +
  facet_grid(~ fld_dim) +
  geom_boxplot() +
  scale_y_continuous(trans = "log")

ggplot(dat, aes(y = KL, x = Fit, fill = Fit)) +
  facet_grid(cols = vars(fld_dim), rows = vars(MODEL)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log")

## REPORTING  RESULTS ##########################################################

plot.res <- 500
baseline_fld_dim <- 100 # shows a typcial results across simulation scenarios
kl.lims <- range(dat$KL[dat$KL < 1e5]) # this better captures some "outliers" identified by default boxplot()
mae.lims <- range(dat$MAE) # this better captures some "outliers" identified by default boxplot()
time.lims <- range(dat$TIME) # likewise exclude these from recorded timing

# subset to just the IDM
plotdat <- dat[dat$MODEL == "IDM" & dat$fld_dim == baseline_fld_dim & dat$Fit != "scampr (opt. k)", ]
plotdat$Fit <- factor(plotdat$Fit, levels = c("INLA", "scampr", "mgcv (def.)", "mgcv (opt.)"))

# for text look at the total computation times
plotdat %>% group_by(Fit) %>% summarise(sum(TIME))

# set the colours to be used in plots - differently if requiring grey-scale or not
make.grey.scale <- FALSE
if (make.grey.scale) {
  fill_cols <- rep("black", 4)
} else {
  fill_cols <- c("darkorange1", "dodgerblue3", "aquamarine4", "darkorchid4")
}

png(filename = paste0(getwd(), "/Figures/baseline_comparison.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
# layout(mat = matrix(c(1:3, 3), nrow = 2, ncol = 2), heights = c(0.5,0.5), widths = c(0.8,0.2))
par(mfrow = c(2, 1))

# accuracy
par(mar = c(1.5, 4.1, 2.1, 0))
boxplot(KL ~ Fit, data = plotdat,
        log = "y", col = alpha(fill_cols, 0.25), xlab = "",
        xaxt = "n", ylab = "", border = fill_cols, yaxt = "n"#, outline = F
)
axis(2, at =c(50, 150, 500, 1500, 5000) , labels = c(50, 150, 500, 1500, 5000))
title(ylab = expression(paste(D[KL], "(", mu, " || ", hat(mu), ")")), cex.lab = 1, line = 3)
mtext("A", side = 3, line = 1.5, adj = 0, padj = 1)

# timing
par(mar = c(2, 4.1, 1.5, 0))
boxplot(TIME ~ Fit, data = plotdat,
        log = "y", col = alpha(fill_cols, 0.25), xaxt = "n", #ylim = time.lims, 
        ylab = "", xlab = "", border = fill_cols, outline = F, yaxt = "n"
)
title(ylab = expression(" Comp. Time"), cex.lab = 1, line = 3)
axis(2, at = c(2, 5, 10, 30, 60) , labels = c("2''", "5''", "10''", "30''", "1'"))
mtext("B", side = 3, line = 1.5, adj = 0, padj = 1)
axis(1, at = 1:length(levels(plotdat$Fit)) , labels = levels(plotdat$Fit))
par(xpd = T)
# WHEN COMP. TIME IS PLOT B
if (!make.grey.scale) {
  points(x = 1:length(levels(plotdat$Fit)) - c(0.25, 0.35, 0.45, 0.45), y = rep(0, length(levels(plotdat$Fit))),
         pch = 22, col = fill_cols, bg = alpha(fill_cols, 0.25))
}
# # WHEN KL Div. IS PLOT B
# if (!make.grey.scale) {
#   points(x = 1:length(levels(plotdat$Fit)) - c(0.25, 0.35, 0.45, 0.45), y = rep(19, length(levels(plotdat$Fit))),
#          pch = 22, col = fill_cols, bg = alpha(fill_cols, 0.25))
# }
# par(mar = rep(0, 4), xpd = F)
# plot(1, type = "n", axes = F) # dummy
# legend("center", title = "IDM fitted via:",
#        legend = levels(plotdat$Fit),
#        col = fill_cols, pch = 15, bty = "n")
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

dev.off()

# comparison over different basis dimensions ###################################

tmp1 <- dat[,c("fld_dim", "sim", "Fit", "MODEL", "KL")]
tmp2 <- dat[,c("fld_dim", "sim", "Fit", "MODEL", "TIME")]
tmp3 <- dat[,c("fld_dim", "sim", "Fit", "MODEL", "MAE")]
colnames(tmp1)[5] <- "value"
colnames(tmp2)[5] <- "value"
colnames(tmp3)[5] <- "value"
tmp1$metric <- "KL"
tmp2$metric <-"Comp. Time (sec)"
tmp3$metric <- "MAE"
rdat <- rbind(tmp1, tmp2, tmp3)

ggplot(dat, aes(y = KL, x = Fit, fill = Fit)) +
  facet_grid(cols = vars(fld_dim), rows = vars(MODEL)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log")

ggplot(rdat %>% filter(MODEL == "IDM"), aes(y = value, x = Fit, fill = Fit)) +
  facet_grid(cols = vars(fld_dim), rows = vars(metric), scale = "free_y") +
  geom_boxplot() +
  scale_y_continuous(trans = "log")

ks <- sort(unique(dat$fld_dim))
idat <- dat %>% filter(MODEL == "IDM" & is.finite(KL) & Fit != "scampr (opt. k)")
idat$Fit <- factor(idat$Fit, levels = c("INLA", "scampr", "mgcv (def.)", "mgcv (opt.)"))
beta.est <- idat %>% group_by(fld_dim, Fit) %>% summarise(rmse = sqrt(mean((BETA_ENV - 1)^2)), std = sd((BETA_ENV - 1)^2))
beta.est$upper <- beta.est$rmse + beta.est$std
beta.est$lower <- beta.est$rmse - beta.est$std
kl.lims <- range(idat$KL)
time.lims <- range(idat$TIME)
beta.lims <- range(beta.est$upper, beta.est$lower)
left_base_mar <- 4.1
left_mars <- c(left_base_mar, left_base_mar * (2/3), left_base_mar / 3, 0)
right_mars <- rev(left_mars)
fit_cols <- c("darkorange1", "dodgerblue3", "aquamarine4", "darkorchid4")

png(filename = paste0(getwd(), "/Figures/comparison_of_basis_dim.png"), width = 6.2 * plot.res, height = 5.2 * plot.res, res = plot.res)
layout(mat = rbind(matrix(1:12, nrow = 3, ncol = 4), rep(13, 4)), heights = c(rep(0.9/3, 3), 0.1), widths = rep(1/4, 4))
for (k in 1:length(ks)) {
  # timing
  par(mar = c(0, left_mars[k], 2.1, right_mars[k]))
  boxplot(TIME ~ Fit, data = idat[idat$fld_dim == ks[k], ],
          log = "y",
          col = alpha(fit_cols, 0.15), xaxt = "n", ylim = time.lims,
          yaxt = "n", ylab = "", border = fit_cols, outline = F
  )
  if (k == 1) {axis(2, at = c(5,20,60,180,480) , labels = c("5''","20''","1'","3'", "8'"))}
  if (k == 1) {title(ylab = expression(paste("Comp. Time")), cex.lab = 1.25, line = 2.5)}
  mtext(paste("k = ", ks[k]), side = 3, line = 1.5, adj = 0, padj = 1)
  # accuracy KLD
  par(mar = c(0, left_mars[k], 2.1, right_mars[k]))
  boxplot(KL ~ Fit, data = idat[idat$fld_dim == ks[k], ],
          log = "y", col = alpha(fit_cols, 0.15), xlab = "", ylim = kl.lims,
          xaxt = "n", yaxt = "n", ylab = "", border = fit_cols
  )
  if (k == 1) {axis(2, at = c(50, 150, 500, 1500) , labels = c(50, 150, 500, 1500))}
  if (k == 1) {title(ylab = expression(paste(D[KL], "(", mu, " || ", hat(mu), ")")), cex.lab = 1.25, line = 2.5)}
  # beta estimation
  par(mar = c(0, left_mars[k], 2.1, right_mars[k]))
  plot(1:length(beta.est$Fit[beta.est$fld_dim == ks[k]]), beta.est$rmse[beta.est$fld_dim == ks[k]],
       ylim = beta.lims, pch = 16, col = fit_cols, yaxt = "n", xaxt = "n", ylab = "", xlab = "",
       xlim = c(0.5, length(beta.est$Fit[beta.est$fld_dim == ks[k]]) + 0.5))
  arrows(x0 = 1:length(beta.est$Fit[beta.est$fld_dim == ks[k]]),
         y0 = beta.est$upper[beta.est$fld_dim == ks[k]],
         y1 = beta.est$lower[beta.est$fld_dim == ks[k]],
         angle = 90, length = 0.05, code = 3, col = fit_cols)
  if (k == 1) {axis(2, at = seq(0.1, 0.3, 0.1) , labels = seq(0.1, 0.3, 0.1))}
  if (k == 1) {title(ylab = expression(paste("RMSE ", hat(beta)[1])), cex.lab = 1.25, line = 2.5)}
}
par(mar = c(0,0,1,0))
plot(1, type = "n", axes = F) # dummy
legend("center", title = "ISDM fitted via:", legend = levels(idat$Fit), horiz = T,
       col = fit_cols, pch = 22, pt.bg = alpha(fit_cols, 0.25), bty = "n", xpd = T, cex = 1.35)
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

dev.off()

## Results for short range latent effects to show when mgcv opt. is necessary ##

res_type <- "_all_short"

# initialise the result storage
dat_short <- NULL

res.list <- list()
res.objs <- list.files(paste0(home.wd, "/Results", res_type))[grepl("res_", list.files(paste0(home.wd, "/Results", res_type)), fixed = T)]
for (job in 1:length(res.objs)) {
  load(paste0(home.wd, "/Results", res_type, "/", res.objs[job]))
  res.list[[job]] <- res_tab
}
dat_short <- do.call(rbind, res.list)
rm(res.list)

# create basis function identifier
dat_short$BASIS_TYPE <- "Optimised"
dat_short$BASIS_TYPE[dat_short$fit_model %in% c("MGCV FIXED", "SCAMPR FIXED")] <- "Fixed"

# adjust the fit to reflect only the software
dat_short$FIT[dat_short$FIT == "SCAMPR FIXED"] <- "SCAMPR"
# dat$Model <- factor(paste0(dat$FIT, " ", dat$MODEL),
#                     levels = c("INLA PA", "SCAMPR PA", "SCAMPR2 PA", "INLA IDM", "SCAMPR IDM", "SCAMPR2 IDM", "INLA PO", "SCAMPR PO", "SCAMPR2 PO"))
dat_short$fit_model[dat_short$fit_model %in% c("INLA PA", "INLA PO", "INLA IDM")] <- "INLA"
dat_short$Fit <- factor(dat_short$fit_model,
                        levels = c("INLA", "SCAMPR FIXED", "SCAMPR", "MGCV FIXED", "MGCV"),
                        labels = c("INLA", "scampr", "scampr (opt. k)","mgcv (def.)", "mgcv (opt.)")
)

# get the rates of failure to converge or other problems
fail.rates.timing_short <- dat_short %>% group_by(MODEL, FIT, BASIS_TYPE) %>%
  summarise(poor.conv.kl = sum(!(KL < 1e5 & !is.na(KL))), poor.conv.mae = sum(!(MAE < 1e5 & !is.na(MAE))), cpu = mean(TIME), fails = 100 - length(unique(sim)))

idat_short <- dat_short %>% filter(MODEL == "IDM" & is.finite(KL) & Fit != "scampr (opt. k)")
idat_short$Fit <- factor(idat_short$Fit, levels = c("INLA", "scampr", "mgcv (def.)", "mgcv (opt.)"))
kl.lims_short <- range(idat_short$KL)

png(filename = paste0(getwd(), "/Figures/comparison_of_long_and_short_range.png"), width = 6.2 * plot.res, height = 5.2 * plot.res, res = plot.res)
layout(mat = rbind(matrix(1:8, nrow = 2, ncol = 4), rep(9, 4)), heights = c(rep(0.9/2, 2), 0.1), widths = rep(1/4, 4))
for (k in 1:length(ks)) {
  # accuracy KLD - long range latent
  par(mar = c(0, left_mars[k], 2.1, right_mars[k]))
  boxplot(KL ~ Fit, data = idat[idat$fld_dim == ks[k], ],
          log = "y", col = alpha(fit_cols, 0.15), xlab = "", ylim = kl.lims,
          xaxt = "n", yaxt = "n", ylab = "", border = fit_cols
  )
  if (k == 1) {axis(2, at = c(50, 150, 500, 1500) , labels = c(50, 150, 500, 1500))}
  if (k == 1) {title(ylab = expression(paste(D[KL], "(", mu, " || ", hat(mu), ")")), cex.lab = 1.25, line = 2.5)}
  mtext(paste("k = ", ks[k]), side = 3, line = 1.5, adj = 0, padj = 1)
  if (k == 4) {mtext("Coarse scale latent range", side = 4, line = 0, adj = Inf, padj = 1)}
  if (k == 4) {mtext(expression(rho[X]<rho[xi]), side = 4, line = 2, adj = Inf, padj = 1)}
  # accuracy KLD - short range latent
  par(mar = c(0, left_mars[k], 2.1, right_mars[k]))
  boxplot(KL ~ Fit, data = idat_short[idat_short$fld_dim == ks[k], ],
          log = "y", col = alpha(fit_cols, 0.15), xlab = "", ylim = kl.lims_short,
          xaxt = "n", yaxt = "n", ylab = "", border = fit_cols
  )
  if (k == 1) {axis(2, at = c(50, 150, 500, 1500) , labels = c(50, 150, 500, 1500))}
  if (k == 1) {title(ylab = expression(paste(D[KL], "(", mu, " || ", hat(mu), ")")), cex.lab = 1.25, line = 2.5)}
  if (k == 4) {mtext("Fine scale latent range", side = 4, line = 0, adj = Inf, padj = 1)}
  if (k == 4) {mtext(expression(rho[X]>rho[xi]), side = 4, line = 2, adj = Inf, padj = 1)}
}
par(mar = c(0,0,1,0))
plot(1, type = "n", axes = F) # dummy
legend("center", title = "IDM fitted via:", legend = levels(idat_short$Fit), horiz = T,
       col = fit_cols, pch = 22, pt.bg = alpha(fit_cols, 0.25), bty = "n", xpd = T, cex = 1.35)
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

dev.off()

# for reporting in text and appendices
plotdat[plotdat$KL > 1e5, ]
plotdat[plotdat$TIME > time.lims[2] | plotdat$TIME < time.lims[1], ]
sim1.tab <- dat %>% group_by(FIT, MODEL, scenario) %>%
  summarise(`kld` = round(median(KL, na.rm = T), 2), `Comp. Time` =round(mean(ALL_TIME), 2), `# kld > 100,000` = sum(KL > 1e5 | is.na(KL)), `# Fit Fail` = as.integer(100 - length(unique(sim))))
library(xtable)
print(
  xtable(
    sim1.tab, label = "append:tab:comp_sim", digits = 2, caption.placement = "top", caption = "Result comparisons for models fit to the PA, PO data separately and jointly (PA, PO and IDM resp.) using either texttt{INLA} or our proposed approach (texttt{scampr} using either optimised basis functions, A, or a set default configuration of $400$, B) over $100$ simulations. Scenarios include both approaches using either a regular grid of domain points, or the default INLA mesh points, to approximate the spatial integral within the PO likelihood component ($n_{text{quad}}=$ 10,000 or 254 resp.). Column emph{$kld$} shows the median Kullback-Leibler divergence from the fitted to true mean abundance rate (where values closer to zero means a more accurate model fit). Column emph{Comp. Time (with basis opt.)} describes median seconds taken to fit the models (the number in brackets additionally includes the seconds taken to optimise the basis function configuration). Column emph{# $kld > $ 100,000} describes the number of simulations for which each model fitted the true mean abundance rate with a Kullback-Leibler divergence of 100,000 or more --- this indicates poor convergence of the model fitting routine. Column emph{# Fit Fail} describes the number of times each model failed to converge."
  ) , include.rownames = F
)