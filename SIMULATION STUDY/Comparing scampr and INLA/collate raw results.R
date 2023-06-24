## Analyse Scenario ##

library(dplyr)
library(ggplot2)

home.wd <- getwd()
# Get the job array
tab <- read.csv("job_array.csv")

# initialise the result storage
dat <- NULL
dat0 <- NULL

res.folders <- c("Results")

for (which.results in 1:length(res.folders)) {
  # initialise the result storage
  dat <- NULL
  dat0 <- NULL
  
  # change into appropriate result folder
  setwd(paste(home.wd, res.folders[which.results], sep = "/"))
  # get a list of all the individual result files
  res.list <- list.files()[grepl("res_", list.files(), fixed = T)]
  # inner loop through individual sim-by-model files
  for (job in res.list) {
    # load in the individual simulation results
    load(job)
    # add to the data
    dat <- rbind(dat, res_tab)
    rm(res_tab)
  }
}
setwd(home.wd)

# set up a factor describing the scenario
dat$scenario <- factor(dat$quad,
                       levels = c("full", "mesh"),
                       labels = c("Entire Grid: $n_{quad}$ = 10,000", "INLA Mesh: $n_{quad}$ = 254"))
# create a new model factor that describes the model and method of fit
dat$FIT[grepl("FIXED", dat$fit_model, fixed = T)] <- paste0(dat$FIT[grepl("FIXED", dat$fit_model, fixed = T)], "2")
dat$Model <- factor(paste0(dat$FIT, " ", dat$MODEL),
                    levels = c("INLA PA", "SCAMPR PA", "SCAMPR2 PA", "INLA IDM", "SCAMPR IDM", "SCAMPR2 IDM", "INLA PO", "SCAMPR PO", "SCAMPR2 PO"))

# get the rates of failure to converge
fail.rates.timing <- dat %>% group_by(Model, scenario) %>%
  summarise(poor.conv.kl = sum(!(KL < 1e5 & !is.na(KL))), poor.conv.mae = sum(!(MAE < 1e5 & !is.na(MAE))), cpu = mean(TIME), fails = 100 - length(unique(sim)))
# dat %>% group_by(Model, scenario) %>%
#   summarise(sim_missing = (1:100)[!(1:100) %in% sim])

# remove failed fits for plotting
plotdat <- dat[dat$Model %in% c("INLA IDM", "SCAMPR IDM", "SCAMPR2 IDM"), ]
plotdat$Model <- factor(plotdat$Model,
                        levels = c("INLA IDM", "SCAMPR2 IDM", "SCAMPR IDM")
)
plotdat$Fit <- factor(plotdat$Model,
                        levels = c("INLA IDM", "SCAMPR2 IDM", "SCAMPR IDM"),
                        labels = c("INLA", "scampr (fix)", "scampr (opt)")
)
ploterr <- fail.rates.timing[fail.rates.timing$Model %in% c("INLA IDM", "SCAMPR IDM", "SCAMPR2 IDM"), ]
ploterr$Model <- factor(ploterr$Model)
ploterr$fail.rate <- ploterr$fails / 100

# set the plotting parameters

fill_cols <- c("darkorange1","dodgerblue4","dodgerblue1")
kl.lims <- range(plotdat$KL[plotdat$KL < 10000]) # this better captures some "outliers" identified by default boxplot()
time.lims <- range(plotdat$TIME[plotdat$KL < 10000]) # likewise exclude these from recorded timing
err.lims <- range(ploterr$fails)# / 100
col1st.width <- 0.35
collast.width <- 0.10
plot.res <- 500

png(filename = paste0(getwd(), "/inla_v_scampr.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
layout(mat = matrix(c(1:5, 3, 6:7, 3), nrow = 3, ncol = length(levels(plotdat$scenario)) + 1, byrow = TRUE), heights = c(0.4,0.4,0.2), widths = c(col1st.width, (1-(col1st.width + collast.width))/2, collast.width))

# computation times
for (i in 1:length(levels(plotdat$scenario))) {
  if (i == 1) {
    par(mar = c(1.5, 4.1, 3, 0))
    boxplot(ALL_TIME ~ Model, data = plotdat[plotdat$scenario == levels(plotdat$scenario)[i], ],
            log = "y", col = fill_cols, xaxt = "n", #ylim = time.lims, 
            ylab = "", border = "black", outline = F
    )
    title(ylab = "Comp. Time (sec)", cex.lab = 1.3, line = 2,
          main = c(expression(bold(n[quad]=="10,000")), expression(bold(n[quad]=="254")))[i]
    )
    mtext("A", side = 3, line = 2, adj = 0, padj = 1)
  } else {
    par(mar = c(1.5, 0, 3, 0))
    boxplot(ALL_TIME ~ Model, data = plotdat[plotdat$scenario == levels(plotdat$scenario)[i], ],
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