## Analyse Scenario ##

library(dplyr)
library(ggplot2)

# toggle for all simulations or just those for the main paper
res_type <- ""
# res_type <- "_all"

home.wd <- getwd()
# Get the job array
tab <- read.csv(paste0("job_array", res_type, ".csv"))

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
                      labels = c("INLA", "scampr (fix)", "scampr (opt)","mgcv (fix)", "mgcv (opt)")
)


# get the rates of failure to converge or other problems
fail.rates.timing <- dat %>% group_by(MODEL, FIT, BASIS_TYPE) %>%
  summarise(poor.conv.kl = sum(!(KL < 1e5 & !is.na(KL))), poor.conv.mae = sum(!(MAE < 1e5 & !is.na(MAE))), cpu = mean(TIME), fails = 100 - length(unique(sim)))
# ploterr <- fail.rates.timing[fail.rates.timing$Model %in% c("INLA IDM", "SCAMPR IDM", "SCAMPR2 IDM"), ]
# ploterr$Model <- factor(ploterr$Model)
# ploterr$fail.rate <- ploterr$fails / 100

ggplot(dat, aes(y = KL, x = Fit)) +
  facet_wrap(~ MODEL + BASIS_TYPE, scales = "free") +
  geom_boxplot()

# subset to just the IDM
plotdat <- dat[dat$MODEL == "IDM", ]

# for text look at the total computation times
plotdat %>% group_by(Fit) %>% summarise(sum(TIME))

# set the plotting parameters

fill_cols <- c("darkorange1","dodgerblue4","dodgerblue1", "aquamarine4", "aquamarine2")
kl.lims <- range(dat$KL) # this better captures some "outliers" identified by default boxplot()
mae.lims <- range(dat$MAE) # this better captures some "outliers" identified by default boxplot()
time.lims <- range(dat$TIME) # likewise exclude these from recorded timing
# err.lims <- range(ploterr$fails)# / 100
plot.res <- 500

png(filename = paste0(getwd(), "/Figures/baseline_comparison.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
layout(mat = matrix(c(1:3, 3), nrow = 2, ncol = 2), heights = c(0.5,0.5), widths = c(0.8,0.2))

# timing
par(mar = c(1.5, 3.1, 2.1, 0))
boxplot(TIME ~ Fit, data = plotdat,
        log = "y", col = fill_cols, xaxt = "n", #ylim = time.lims, 
        ylab = "", border = "black", outline = F
)
title(ylab = "Avg. Comp. Time (sec)", cex.lab = 1, line = 2)
      # main = c(expression(bold(n[quad]=="10,000")), expression(bold(n[quad]=="254")))[i])
mtext("A", side = 3, line = 2, adj = 0, padj = 1)

# accuracy
par(mar = c(1.5, 3.1, 1.5, 0))
boxplot(KL ~ Fit, data = plotdat,
        log = "y", col = fill_cols, xlab = "",
        xaxt = "n", yaxt = "n", ylab = "", border = "black"#, outline = F
)
axis(2, at =c(50, 150, 500, 1500, 5000) , labels = c(50, 150, 500, 1500, 5000))
title(ylab = expression(paste(D[KL], "(", mu, " || ", hat(mu), ")")), cex.lab = 1, line = 2)
mtext("B", side = 3, line = 2, adj = 0, padj = 1)

par(mar = rep(0, 4))
plot(1, type = "n", axes = F) # dummy
legend("center", title = "IDM fitted via:",
       legend = levels(plotdat$Fit),
       col = fill_cols, pch = 15, bty = "n")
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))

dev.off()

## Plot for all models in the baseline scenario ################################

png(filename = paste0(getwd(), "/Figures/baseline_comparison_of_models.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
layout(mat = matrix(c(1:9, rep(10, 3)), nrow = 3, ncol = 4), heights = rep(1/3, 3), widths = c(rep(0.85/3, 3),0.15))

for (i in 1:3) {
  
  # set some panel-specific stuff
  tmp.mod <- c("IDM", "PA", "PO")[i]
  tmp.mar <- list(c(3.1, 0), c(1.55, 1.55), c(0, 3.1))
  
  # timing
  par(mar = c(0, tmp.mar[[i]][1], 2.1, tmp.mar[[i]][2]))
  boxplot(TIME ~ Fit, data = dat[dat$MODEL == tmp.mod, ],
          log = "y", col = fill_cols, xaxt = "n", ylim = time.lims,
          yaxt = "n", ylab = "", border = "black", outline = F
  )
  if (i == 1) {axis(2, at = c(0.2, 0.5, 2, 5, 20) , labels = c("0.2", "0.5", "2", "5", "20"))}
  title(ylab = "Avg. Comp. Time (sec)", cex.lab = 1, line = 2)
  mtext(tmp.mod, side = 3, line = 1.5, adj = 0, padj = 1)
  
  # accuracy KLD
  par(mar = c(0, tmp.mar[[i]][1], 2.1, tmp.mar[[i]][2]))
  boxplot(KL ~ Fit, data = dat[dat$MODEL == tmp.mod, ],
          log = "y", col = fill_cols, xlab = "", ylim = kl.lims,
          xaxt = "n", yaxt = "n", ylab = "", border = "black"#, outline = F
  )
  if (i == 1) {axis(2, at = c(50, 150, 500, 2500, 10000) , labels = c(50, 150, 500, 2500, "10K"))}
  title(ylab = expression(paste(D[KL], "(", mu, " || ", hat(mu), ")")), cex.lab = 1, line = 2)
  
  # accuracy MAE
  par(mar = c(0, tmp.mar[[i]][1], 2.1, tmp.mar[[i]][2]))
  boxplot(MAE ~ Fit, data = dat[dat$MODEL == tmp.mod, ], log = "y", ylim = mae.lims, 
          col = fill_cols, xlab = "", xaxt = "n", ylab = "", border = "black", yaxt = "n"
  )
  if (i == 1) {axis(2, at = c(0.05, 0.15, 0.4, 1), labels = c(0.05, 0.15, 0.4, 1))}
  title(ylab = expression(paste(MAE, " (", mu, " - ", hat(mu), ")")), cex.lab = 1, line = 2)
  

}
par(mar = rep(0, 4))
plot(1, type = "n", axes = F) # dummy
legend("center", title = "Fitted via:",
       legend = levels(plotdat$Fit),
       col = fill_cols, pch = 15, bty = "n", xpd = T, cex = 1.25)
par(mfrow = c(1,1), mar = c(5.1,4.1,4.1,2.1))


dev.off()

# computation times
for (i in 1:length(levels(plotdat$scenario))) {
  if (i == 1) {
    par(mar = c(1.5, 4.1, 3, 0))
    boxplot(TIME ~ Model, data = plotdat[plotdat$scenario == levels(plotdat$scenario)[i], ],
            log = "y", col = fill_cols, xaxt = "n", #ylim = time.lims, 
            ylab = "", border = "black", outline = F
    )
    title(ylab = "Comp. Time (sec)", cex.lab = 1.3, line = 2,
          main = c(expression(bold(n[quad]=="10,000")), expression(bold(n[quad]=="254")))[i]
    )
    mtext("A", side = 3, line = 2, adj = 0, padj = 1)
  } else {
    par(mar = c(1.5, 0, 3, 0))
    boxplot(TIME ~ Model, data = plotdat[plotdat$scenario == levels(plotdat$scenario)[i], ],
            log = "y", col = fill_cols, xaxt = "n", #ylim = time.lims, 
            yaxt = "n", border = "black", outline = F
    )
    title(ylab = "Comp. Time (sec)", cex.lab = 1.3, line = 2,
          main = c(expression(bold(n[quad]=="10,000")), expression(bold(n[quad]=="254")))[i]
    )
  }
}
par(mar = c(0, 0, 0, 0))

# legend
plot(1, type = "n", axes = F) # dummy
legend("center", title = "IDM fitted via:",
       legend = levels(plotdat$Fit),
       col = fill_cols, pch = 15, bty = "n")

# accuracy
for (i in 1:length(levels(plotdat$scenario))) {
  if (i == 1) {
    par(mar = c(1.5, 4.1, 1.5, 0))
    boxplot(KL ~ Model, data = plotdat[plotdat$scenario == levels(plotdat$scenario)[i], ],
            log = "y", col = fill_cols, ylim = kl.lims, xlab = "",
            xaxt = "n", yaxt = "n", ylab = "", border = "black"#, outline = F
    )
    axis(2, at =c(50, 150, 500, 1500, 5000) , labels = c(50, 150, 500, 1500, 5000))
    title(ylab = expression(paste(D[KL], "(", mu, " || ", hat(mu), ")")), cex.lab = 1.3, line = 2)
    mtext("B", side = 3, line = 2, adj = 0, padj = 1)
  } else {
    par(mar = c(1.5, 0, 1.5, 0))
    boxplot(KL ~ Model, data = plotdat[plotdat$scenario == levels(plotdat$scenario)[i], ],
            log = "y", col = fill_cols, main = "", xlab = "", ylim = kl.lims, 
            xaxt = "n", yaxt = "n", border = "black"#, outline = F
    )
  }
}


png(filename = paste0(getwd(), "/inla_v_scampr.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
layout(mat = matrix(c(1:5, 3, 6:7, 3), nrow = 3, ncol = length(levels(plotdat$scenario)) + 1, byrow = TRUE), heights = c(0.4,0.4,0.2), widths = c(col1st.width, (1-(col1st.width + collast.width))/2, collast.width))

# computation times
for (i in 1:length(levels(plotdat$scenario))) {
  if (i == 1) {
    par(mar = c(1.5, 4.1, 3, 0))
    boxplot(TIME ~ Model, data = plotdat[plotdat$scenario == levels(plotdat$scenario)[i], ],
            log = "y", col = fill_cols, xaxt = "n", #ylim = time.lims, 
            ylab = "", border = "black", outline = F
    )
    title(ylab = "Comp. Time (sec)", cex.lab = 1.3, line = 2,
          main = c(expression(bold(n[quad]=="10,000")), expression(bold(n[quad]=="254")))[i]
    )
    mtext("A", side = 3, line = 2, adj = 0, padj = 1)
  } else {
    par(mar = c(1.5, 0, 3, 0))
    boxplot(TIME ~ Model, data = plotdat[plotdat$scenario == levels(plotdat$scenario)[i], ],
            log = "y", col = fill_cols, xaxt = "n", #ylim = time.lims, 
            yaxt = "n", border = "black", outline = F
    )
    title(ylab = "Comp. Time (sec)", cex.lab = 1.3, line = 2,
          main = c(expression(bold(n[quad]=="10,000")), expression(bold(n[quad]=="254")))[i]
    )
  }
}
par(mar = c(0, 0, 0, 0))

# legend
plot(1, type = "n", axes = F) # dummy
legend("center", title = "IDM fitted via:",
       legend = levels(plotdat$Fit),
       col = fill_cols, pch = 15, bty = "n")

# accuracy
for (i in 1:length(levels(plotdat$scenario))) {
  if (i == 1) {
    par(mar = c(1.5, 4.1, 1.5, 0))
    boxplot(KL ~ Model, data = plotdat[plotdat$scenario == levels(plotdat$scenario)[i], ],
            log = "y", col = fill_cols, ylim = kl.lims, xlab = "",
            xaxt = "n", yaxt = "n", ylab = "", border = "black"#, outline = F
    )
    axis(2, at =c(50, 150, 500, 1500, 5000) , labels = c(50, 150, 500, 1500, 5000))
    title(ylab = expression(paste(D[KL], "(", mu, " || ", hat(mu), ")")), cex.lab = 1.3, line = 2)
    mtext("B", side = 3, line = 2, adj = 0, padj = 1)
  } else {
    par(mar = c(1.5, 0, 1.5, 0))
    boxplot(KL ~ Model, data = plotdat[plotdat$scenario == levels(plotdat$scenario)[i], ],
            log = "y", col = fill_cols, main = "", xlab = "", ylim = kl.lims, 
            xaxt = "n", yaxt = "n", border = "black"#, outline = F
    )
  }
}

# error rates
for (i in 1:length(levels(ploterr$scenario))) {
  if (i == 1) {
    par(mar = c(2.1, 4.1, 1.5, 0))
    barplot(fails ~ Model, data =  ploterr[ploterr$scenario == levels(ploterr$scenario)[i], ],
            col = fill_cols, xaxt = "n", ylab = "", ylim = err.lims
    )
    title(ylab = "% Failed", cex.lab = 1.3, line = 2)
    mtext("C", side = 3, line = 2, adj = 0, padj = 1)
  } else {
    par(mar = c(2.1, 0, 1.5, 0))
    barplot(fails ~ Model, data =  ploterr[ploterr$scenario == levels(ploterr$scenario)[i], ],
            col = fill_cols, xaxt = "n", ylab = "", yaxt = "n", ylim = err.lims
    )
  }
}
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