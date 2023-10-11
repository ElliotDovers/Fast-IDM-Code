## Analyse Scenario ##

home.wd <- getwd()

# initialise the result storage
dat <- NULL

# change into appropriate result folder
setwd(paste(home.wd, "raw", sep = "/"))
# get a list of all the individual result files
res.list <- list.files()[grepl("res_", list.files(), fixed = T)]
# inner loop through individual sim-by-model files
for (job in res.list) {
  # get the job number
  job.no <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(job, ".", fixed = T), function(x){x[1]})), "_", fixed = T), function(x){x[length(x)]})))
  # load in the individual simulation results
  load(job)
  # add in the job number
  res_tab$job <- job.no
  dat <- rbind(dat, res_tab)
  rm(res_tab)
}
setwd(home.wd)

library(dplyr)
library(ggplot2)
library(gridExtra)

# Check out the number of flagged convergences - this are fit failures for scampr
fit.fails_scampr <- dat %>% group_by(Fit, Model) %>% summarise(OK = sum(fit_flag == 0), FAIL = sum(fit_flag == 1))
# Check out the missing simulations for INLA, i.e. fit failures
fit.fails_inla <- dat %>% group_by(Fit, Model) %>% summarise(OK = sum(fit_flag == 0), FAIL = sum(fit_flag == 1))
dat <- dat[dat$fit_flag != 1, ]

# create factors from the model type:
dat$Model <- factor(dat$Model,
                    levels = c('structured', 'unstructured', 'joint', 'unstructuredcov', 'jointcov', 'jointtwo'),
                    labels = c('PA only (A)', 'PO only (B)', 'IDM (C)', 'PO with \nbias \ncovariate (D)', 'IDM with \nbias \ncovariate (E)', 'IDM with \nsecond latent \nfield (F)')
)


# Separate out the models for different comparisons:
dat0 <- dat[dat$Fit %in% c("INLA", "scampr"), ]
dat1 <- dat[!dat$Fit %in% c("INLA", "scampr"), ]
dat0 <- dat[dat$Fit == "scampr", ]
newdat <- data.frame(score = c(dat0$MAE, dat0$KL_div), metric = rep(c("MAE", "D[KL](mu||hat(mu))"), each = nrow(dat0)), model = c(dat0$Model, dat0$Model))
ggplot(data = newdat, aes(y = score)) + geom_boxplot(outlier.shape = NA) +
  facet_grid(rows = vars(metric), cols = vars(model), scales = "free")#, labeller = label_bquote(rows = k~.(as.character(k_fact)), cols = m==.(m))) +

# re-level factors
dat0$Fit <- factor(as.character(dat0$Fit))
dat1$Model <- factor(as.character(dat1$Model))
# dat1$Fit <- factor(dat1$Fit, levels = c('scampr_ll_k_dense_k_bias_coarse', 'scampr_ll_k_coarse_k_bias_dense'),
#                    labels = c("config. 1", "config. 2"))

# calculate the y axis limits for plotting each metric using a boxplot without the outliers:
MAE.ylims <- dat0 %>%
  group_by(Model, Fit) %>%
  summarise(ymin = boxplot.stats(MAE)$stats[1], ymax = boxplot.stats(MAE)$stats[5]) %>%
  ungroup() %>%
  summarise(ymin = min(ymin), ymax = max(ymax))
KL.ylims <- dat0 %>% group_by(Model, Fit) %>% summarise(ymin = boxplot.stats(KL_div)$stats[1],
                                                         ymax = boxplot.stats(KL_div)$stats[5])
TIME.ylims <- dat0 %>% group_by(Model, Fit) %>% summarise(ymin = boxplot.stats(timing)$stats[1],
                                                        ymax = boxplot.stats(timing)$stats[5])

# # set manual colours FROM SIMMONDS ET AL
# manual_colours <- c("orange", "blue", "grey30", "darkblue",  "grey50", "grey80")

# set the plotting parameters
fill_cols <- c("darkorange1", "dodgerblue1")
line_cols <-  c("darkorange4", "dodgerblue4")
mae.lims <- as.vector(unlist(dat0 %>% #filter(Model %in% levels(dat0$Model)[5:6]) %>%
  group_by(Model, Fit) %>%
  summarise(ymin = boxplot.stats(MAE)$stats[1], ymax = boxplot.stats(MAE)$stats[5]) %>%
  ungroup() %>%
  summarise(ymin = min(ymin), ymax = max(ymax))
))
kl.lims <- as.vector(unlist(dat0 %>% #filter(Model %in% levels(dat0$Model)[5:6]) %>%
  group_by(Model, Fit) %>%
  summarise(ymin = boxplot.stats(KL_div)$stats[1], ymax = boxplot.stats(KL_div)$stats[5]) %>%
  ungroup() %>%
  summarise(ymin = min(ymin), ymax = max(ymax))
))
time.lims <- as.vector(unlist(dat0 %>% #filter(Model %in% levels(dat0$Model)[5:6]) %>%
  group_by(Model, Fit) %>%
  summarise(ymin = boxplot.stats(timing)$stats[1], ymax = boxplot.stats(timing)$stats[5]) %>%
  ungroup() %>%
  summarise(ymin = min(ymin), ymax = max(ymax))
))
col1.width <- 0.225
plot.res <- 500



layout(mat = matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE), heights = c(0.375, 0.3, 0.325), widths = c(col1.width, rep((1-col1.width)/3,3)))
par(xpd = T)
for (i in 5:6) {
  if (i == 5) {
    # par(mgp = c(1.75,0.75,0), xpd = F)
    par(mar = c(0, 3.5, 4.1, 0))
    boxplot(MAE ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
            log = "y", col = fill_cols, ylim = mae.lims, main = levels(dat0$Model)[i],
            xaxt = "n", ylab = "", border = line_cols, outline = F
    )
    title(ylab = "MAE", cex.lab = 1.3, line = 2)
  } else {
    par(mar = c(0, 0, 4.1, 0))
    boxplot(MAE ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
            log = "y", col = fill_cols, ylim = mae.lims, main = levels(dat0$Model)[i],
            xaxt = "n", yaxt = "n", border = line_cols, outline = F
    )
  }
}
for (i in 5:6) {
  if (i == 5) {
    # par(mgp = c(1.75,0.75,0), xpd = F)
    par(mar = c(0, 3.5, 0, 0))
    boxplot(KL_div ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
            log = "y", col = fill_cols, ylim = kl.lims, yaxt = "n",
            xaxt = "n", ylab = "", border = line_cols, outline = F
    )
    axis(2, at = c(50, 100, 500, 2000), labels = c(50, 100, 500, 2000))
    title(ylab = expression(paste(D[KL],"(", mu, " || ", hat(mu), ")")), cex.lab = 1.3, line = 2)
  } else {
    par(mar = c(0, 0, 0, 0))
    boxplot(KL_div ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
            log = "y", col = fill_cols, ylim = kl.lims,
            xaxt = "n", yaxt = "n", border = line_cols, outline = F
    )
  }
}
for (i in 5:6) {
  if (i == 5) {
    par(mar = c(2.1, 3.5, 0, 0))
    boxplot(timing ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
            log = "y", col = fill_cols, ylim = time.lims,
            xaxt = "n", ylab = "", border = line_cols, outline = F
    )
    title(ylab = "Comp. Time (secs)", cex.lab = 1.3, line = 2)
  } else {
    par(mar = c(2.1, 0, 0, 0))
    boxplot(timing ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
            log = "y", col = fill_cols, ylim = time.lims,
            xaxt = "n", yaxt = "n", border = line_cols, outline = F
    )
    if (i == 5) {
      # Legend
      legend(x = 0.5, y = time.lims[1], legend = c("INLA"),
             col = line_cols[1], pch = 22, cex = 1.5, pt.bg = fill_cols[1],
             bty = "n", xpd = T, horiz = T)
      # x.intersp = c(0.5, 0.5, 0.5), text.width = c(0.25,0.27,0.26))
    } else if (i == 6) {
      legend(x = 0.5, y = time.lims[1], legend = c("scampr"),
             col = line_cols[2], pch = 22, cex = 1.5, pt.bg = fill_cols[2], 
             bty = "n", xpd = T, horiz = T)
      # x.intersp = c(0.5, 0.5, 0.5), text.width = c(0.25,0.27,0.26))
    }
  }
}


# plot for appendix including KL Divergence metric
png(filename = paste0(getwd(), "/sym_sim_results_poor_bias_covar.png"), width = 6.2 * plot.res, height = 4.2 * plot.res, res = plot.res)
layout(mat = matrix(1:12, nrow = 2, ncol = 6, byrow = TRUE), heights = c(0.375, 0.325), widths = c(col1.width, rep((1-col1.width)/5,5)))
par(xpd = T)
for (i in 1:6) {
  if (i == 1) {
    # par(mgp = c(1.75,0.75,0), xpd = F)
    par(mar = c(0, 3.5, 4.1, 0))
    boxplot(MAE ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
            log = "y", col = fill_cols[2], ylim = mae.lims, main = levels(dat0$Model)[i],
            xaxt = "n", ylab = "", border = line_cols[2], outline = F
    )
    title(ylab = "MAE", cex.lab = 1.3, line = 2)
  } else {
    par(mar = c(0, 0, 4.1, 0))
    boxplot(MAE ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
            log = "y", col = fill_cols[2], ylim = mae.lims, main = levels(dat0$Model)[i],
            xaxt = "n", yaxt = "n", border = line_cols[2], outline = F
    )
  }
}
for (i in 1:6) {
  if (i == 1) {
    # par(mgp = c(1.75,0.75,0), xpd = F)
    par(mar = c(0, 3.5, 0, 0))
    boxplot(KL_div ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
            log = "y", col = fill_cols[2], ylim = kl.lims, yaxt = "n",
            xaxt = "n", ylab = "", border = line_cols[2], outline = F
    )
    axis(2, at = c(50, 100, 500, 2000), labels = c(50, 100, 500, 2000))
    title(ylab = expression(paste(D[KL],"(", mu, " || ", hat(mu), ")")), cex.lab = 1.3, line = 2)
  } else {
    par(mar = c(0, 0, 0, 0))
    boxplot(KL_div ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
            log = "y", col = fill_cols[2], ylim = kl.lims,
            xaxt = "n", yaxt = "n", border = line_cols[2], outline = F
    )
  }
}
# for (i in 1:6) {
#   if (i == 1) {
#     par(mar = c(2.1, 3.5, 0, 0))
#     boxplot(timing ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
#             log = "y", col = fill_cols, ylim = time.lims,
#             xaxt = "n", ylab = "", border = line_cols, outline = F
#     )
#     title(ylab = "Comp. Time (secs)", cex.lab = 1.3, line = 2)
#   } else {
#     par(mar = c(2.1, 0, 0, 0))
#     boxplot(timing ~ Fit, data = dat0[dat0$Model == levels(dat0$Model)[i], ],
#             log = "y", col = fill_cols, ylim = time.lims,
#             xaxt = "n", yaxt = "n", border = line_cols, outline = F
#     )
#     if (i == 3) {
#       # Legend
#       legend(x = 0.5, y = time.lims[1], legend = c("INLA"),
#              col = line_cols[1], pch = 22, cex = 1.5, pt.bg = fill_cols[1],
#              bty = "n", xpd = T, horiz = T)
#       # x.intersp = c(0.5, 0.5, 0.5), text.width = c(0.25,0.27,0.26))
#     } else if (i == 4) {
#       legend(x = 0.5, y = time.lims[1], legend = c("scampr"),
#              col = line_cols[2], pch = 22, cex = 1.5, pt.bg = fill_cols[2], 
#              bty = "n", xpd = T, horiz = T)
#       # x.intersp = c(0.5, 0.5, 0.5), text.width = c(0.25,0.27,0.26))
#     }
#   }
# }
dev.off()

# Secondary Latent Field Model Comparison #

dat_ll <- NULL
for (i in 1:500) {
  if (i %in% unique(dat1$sim)) {
    dat_ll <- rbind(dat_ll, cbind(ll_diff = dat1$ll[dat1$sim == i & dat1$Fit == "scampr_ll_k_dense_k_bias_coarse"] - dat1$ll[dat1$sim == i & dat1$Fit == "scampr_ll_k_coarse_k_bias_dense"],
                                            sim = i, dense_coarse = dat1$ll[dat1$sim == i & dat1$Fit == "scampr_ll_k_dense_k_bias_coarse"], coarse_dense = dat1$ll[dat1$sim == i & dat1$Fit == "scampr_ll_k_coarse_k_bias_dense"]
                                  )
    )
  } else {
    dat_ll <- rbind(dat_ll, cbind(ll_diff = NA,
                                            sim = i, dense_coarse = NA, coarse_dense = NA))
  }
}
dat_ll <- data.frame(dat_ll)

# adjust for plotting label
dat1$Fit <- factor(dat1$Fit, levels = c('scampr_ll_k_dense_k_bias_coarse', 'scampr_ll_k_coarse_k_bias_dense'),
                   labels = c("Config. 1", "Config. 2"))

# set the plotting parameters
fill_cols_ll <- c("darkseagreen", "darkslategray4")
line_cols_ll <- c("darkseagreen4", "darkslategray")
col1.width <- 1/2
point_fill_cols <- rep(fill_cols_ll[1], nrow(dat_ll))
point_fill_cols[dat_ll$ll_diff < 0] <- fill_cols_ll[2]
point_line_cols <- rep(line_cols_ll[1], nrow(dat_ll))
point_line_cols[dat_ll$ll_diff < 0] <- line_cols_ll[2]
# point.pch <- rep(1, nrow(ll_dat_wide))
# point.pch[ll_dat_wide$ll_diff < 0] <- 1

plot.res <- 500
png(filename = paste0(getwd(), "/sym_sim_selection.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
layout(mat = matrix(c(1, 1, 2, 3, 4, 4), nrow = 3, ncol = 2, byrow = TRUE), heights = c(0.05, 0.8, 0.15), widths = rep(0.5, 2))
# par(xpd = F)
par(mar = c(0, 0, 0, 0))
plot(1, 1, type = "n", axes = F)
text(1, 1, "Performance", cex = 2)
par(mar = c(5.1, 4.1, 4.1, 2.1))
# plot(-2*ll_diff ~ sim, dat_ll, xlab = "Simulation #", ylab = expression(paste(Delta, " AIC")), col = point_line_cols, pch = 21,
#      main = "A: In-sample", bg = point_fill_cols)
plot(-2*ll_diff ~ sim, dat_ll, xlab = "Simulation #", ylab = expression(paste(Delta, " AIC (Config. 1 - Config. 2)")), pch = 21,
     main = "A: In-sample")
abline(h = 0, col = "black", lty = "dashed")
par(mar = c(5.1, 4.1, 4.1, 2.1))
boxplot(KL_div ~ Fit, data = dat1, main = "B: Out-of-sample",
        log = "y", col = fill_cols_ll, xlab = "Dual Resolution Basis Functions",
        ylab = expression(paste(D[KL],"(", mu, " || ", hat(mu), ")")), border = line_cols_ll
)
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("center", legend = expression("Config. 1: k = 20x20, "~k[bias]~"= 5x5", "Config. 2: k = 5x5, "~k[bias]~"= 20x20"), col = line_cols_ll, title = "scampr IDM basis function setup:", pch = 22, pt.bg = fill_cols_ll, bty = "n", cex = 1.5)
dev.off()