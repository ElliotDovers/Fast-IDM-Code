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
  # add in the scenario
  res_tab$scen <- paste0("PO", res_tab$scen_po, "_PA", res_tab$scen_pa)
  # add to the data
  dat <- rbind(dat, res_tab)
  # ranges <- rbind(ranges, attr(res_tab, "correlation ranges"))
  rm(res_tab)
}

setwd(home.wd)

# save(list = "dat", file = "collated raw results.RDATA")
load("collated raw results.RDATA")

library(dplyr)
library(ggplot2)

all(paste0("PO", dat$scen_po, "_PA", dat$scen_pa) == dat$scen) # check I have assigned scenarios correctly
dat$scenario <- factor(paste0(dat$scen, "_", dat$env_range))

# check out the number of flagged convergences - this are fit failures for scampr
dat %>% group_by(scen, Model, env_range) %>%
  summarise(FLAGGED = sum(fit_flag == 1), INF_KL = sum(is.infinite(KLdiv)), INF_LL = sum(is.infinite(LL)), POOR_CONV = sum(LL > 1e5 | is.nan(LL)))

# create full results table for Appendix
fit.issues <- dat %>% group_by(scen, env_range, latent_range, Model) %>%
  summarise(N_SIMS = length((sim)), INF_KL = sum(is.infinite(KLdiv) & !(LL > 1e5 | is.nan(LL))), POOR_CONV = sum(LL > 1e5 | is.nan(LL)))

# check that the infinite KLs are not concentrated so that there are not enough simulations
sim.check <- dat %>% filter(is.finite(KLdiv)) %>% group_by(Model, scen, env_range, bias_range) %>%
  summarise(unique_sims = length(unique(sim)), sims = length(sim), avg_k = mean(k), med_k = median(k), avg_k_bias = mean(k_bias), med_k_bias = median(k_bias))
min(sim.check$sims) # 84% is the lowest number returned in any one combination of env, lat and bias ranges

# check the number of basis functions selected for the problem fits
dat %>% filter(is.finite(KLdiv) | is.na(KLdiv)) %>% group_by(Model, scen, env_range, bias_range) %>%
  summarise(unique_sims = length(unique(sim)), sims = length(sim), avg_k = mean(k), med_k = median(k), avg_k_bias = mean(k_bias), med_k_bias = median(k_bias))
# the majority of problem KL are due to poor convergences
table(dat[is.infinite(dat$KLdiv) | is.na(dat$KLdiv), "LL"] > 1e5)

# remove entire sims with infinite or NaN KLdiv (SO TO COMPARE MODELS FAIRLY)
tmp <- dat %>% filter(Model != "PO Only") %>% group_by(scen, env_range, bias_range, sim) %>% summarise(check = sum(is.infinite(KLdiv) | is.na(KLdiv)))
dat$id <- paste(dat$scen, dat$env_range, dat$bias_range, dat$sim, sep = ":")
# dat[dat$id == paste(tmp[which.max(tmp$check), c("scen", "env_range", "bias_range", "sim")], collapse = ":"), ] # sense check
rm.ids <- paste(tmp$scen[tmp$check != 0], tmp$env_range[tmp$check != 0], tmp$bias_range[tmp$check != 0], tmp$sim[tmp$check != 0], sep = ":")
dat <- dat[!dat$id %in% rm.ids, ]

# there are also 2 KLdiv results that are effectively Infinite that are throwing off results
dat[dat$KLdiv > 1e10, ]
tmp <- dat %>% filter(Model != "PO Only") %>% group_by(scen, env_range, bias_range, sim) %>% summarise(check = sum(KLdiv > 1e10))
rm.ids <- paste(tmp$scen[tmp$check != 0], tmp$env_range[tmp$check != 0], tmp$bias_range[tmp$check != 0], tmp$sim[tmp$check != 0], sep = ":")
dat <- dat[!dat$id %in% rm.ids, ]

# there are also 3 KLdiv results for the PA Only model that are clearly poor convergence (when looking at the other model's KL for those sims)
dat[dat$KLdiv > 1e7, ]
tmp <- dat %>% filter(Model != "PO Only") %>% group_by(scen, env_range, bias_range, sim) %>% summarise(check = sum(KLdiv > 1e7))
rm.ids <- paste(tmp$scen[tmp$check != 0], tmp$env_range[tmp$check != 0], tmp$bias_range[tmp$check != 0], tmp$sim[tmp$check != 0], sep = ":")
dat <- dat[!dat$id %in% rm.ids, ]

# factor the environment variable range with nice labels
dat$fixed_ranges <- factor(dat$env_range,
                      levels = c("15",
                                 "25"
                      ),
                      labels = c("Fine scale Environ. - Coarse scale Latent",
                                 "Coarse scale Environ. - Fine scale Latent"
                      )
)

# adjust our metric to be relative to the PA Only model
newdat <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(delta_kl = KLdiv[Model == "IDM"] - KLdiv[Model == "PA Only"])
newdat$scen_pa <- sapply(strsplit(newdat$scen, "_PA", fixed = T), function(x){x[2]})
newdat$Model <- paste0("IDM - N PO = ", sapply(strsplit(sapply(strsplit(newdat$scen, "PO", fixed = T), function(x){x[2]}), "_", fixed = T), function(x){x[1]}))
newdat$b_range <- factor(newdat$bias_range)
newdat$Model <- factor(newdat$Model, levels = c("IDM - N PO = 50", "IDM - N PO = 200", "IDM - N PO = 800"))

# set the plot resolution
plot.res <- 500

# set a scale function to round to KL div. axis values to zero decimal places
scaleFUN <- function(x) sprintf("%.0f", x)

# OUR MAIN RESULT PLOT (using KL div. on the predicted abundance rate)

plotdat <- dat[dat$Model != "PO Only" & dat$fixed_ranges == "Fine scale Environ. - Coarse scale Latent", ]

# re-factor the model to incorporate PO scenario into IDM
plotdat$Model[plotdat$Model == "IDM"] <- paste(plotdat$Model[plotdat$Model == "IDM"], plotdat$scen_po[plotdat$Model == "IDM"], sep = " - ")
plotdat$Model <- factor(plotdat$Model,
                        levels = c("IDM - 50", "IDM - 200", "IDM - 800", "PA Only"),
                        labels =  c("IDM - N PO = 50", "IDM - N PO = 200", "IDM - N PO = 800", "PA only")
)


model_cols <- c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4", "darkblue")

png(filename = paste0(home.wd, "/Figures/our_sim_results_smoothers.png"), width = 6.3 * plot.res, height = 3.5 * plot.res, res = plot.res)
ggplot(data = plotdat, aes(x = bias_range, y = KLdiv, col = Model)) +
  geom_smooth(aes(fill = Model)) + facet_wrap(~ scen_pa, scales = "free_y", labeller = label_bquote(E*"["*n[PA]*"]"==.(as.numeric(as.character(scen_pa))))) +
  scale_x_continuous(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s)))) +
  scale_y_continuous(name = expression(bold(paste(D[KL],"(", mu, " || ", hat(mu), ")"))), breaks = c(seq(55, 80, by = 5), seq(140, 190, by = 10), seq(4e3, 10e3, by = 2e3), 15e3, 20e3), trans = "log", labels=scaleFUN) +
  scale_color_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  scale_fill_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10),
        legend.text=element_text(size=10), legend.text.align = 0, plot.title = element_text(hjust = 0.5))
dev.off()

## Boxplots relative to PA Only ##

# set the ranges for each n[PA]
tmp <- boxplot(delta_kl ~ scen_pa, data = newdat)
ranges <- data.frame(apply(tmp$stats, 2, function(x){c(floor(x[1]), ceiling(x[5]))}))
colnames(ranges) <- c("50", "200", "800")
newdat$scen_pa <- factor(newdat$scen_pa, levels = c(50, 200, 800))

mod_cols <- c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4")

p1 <- ggplot(data = newdat %>% filter(scen_pa == 50 & b_range == "19"), aes(x = Model, y = delta_kl, color = Model, fill = Model)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected Amount of Presence-Only Data Integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_y_continuous(
    # limits = ranges$`200`,
    name = expression(bold(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")"))),
    # trans = "log",
    # breaks = c(20, 50, 150),
    labels=scaleFUN) +
  coord_cartesian(ylim = c(-100, 50)) +
  scale_color_manual(values = mod_cols) +
  scale_fill_manual(values = alpha(mod_cols, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = 6),
        plot.margin = margin(t = 20,  # Top margin
                             r = 15,  # Right margin
                             b = 12,  # Bottom margin
                             l = 11)  # Left margin
  ) +
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 50"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p2 <- ggplot(data = newdat %>% filter(scen_pa == 200 & b_range == "19"), aes(x = Model, y = delta_kl, color = Model, fill = Model)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected Amount of Presence-Only Data Integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_y_continuous(
    # limits = ranges$`200`,
    name = expression(bold(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")"))),
    # trans = "log",
    # breaks = c(20, 50, 150),
    labels=scaleFUN) +
  coord_cartesian(ylim = c(-180, 90)) +
  scale_color_manual(values = mod_cols) +
  scale_fill_manual(values = alpha(mod_cols, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = 6),
        plot.margin = margin(t = 12,  # Top margin
                             r = 15,  # Right margin
                             b = 20,  # Bottom margin
                             l = 11)  # Left margin
  ) + 
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 200"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p3 <- ggplot(data = newdat %>% filter(scen_pa == 800 & b_range == "19"), aes(x = Model, y = delta_kl, col = Model, fill = Model)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected Amount of Presence-Only Data Integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_y_continuous(
    # limits = ranges$`800`,
    name = expression(bold(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")"))),
    # trans = "log",
    # breaks = c(20, 50, 150),
    labels=scaleFUN) +
  coord_cartesian(ylim = c(-26000, 14000)) +
  scale_color_manual(values = mod_cols) +
  scale_fill_manual(values = alpha(mod_cols, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10), axis.title.x = element_text(margin = margin(t = 15)),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none",
        plot.margin = margin(t = 0,  # Top margin
                             r = 15,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 800"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

library(gridExtra)
png(filename = paste0(home.wd, "/Figures/our_sim_results_boxplots.png"), width = 4.2 * plot.res, height = 6.5 * plot.res, res = plot.res)
grid.arrange(p1,p2,p3)
dev.off()

### Appendix plot - all results

p1 <- ggplot(data = newdat %>% filter(scen_pa == 50), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s)))) +
  scale_y_continuous(
    # limits = ranges$`200`,
    name = expression(bold(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")"))),
    # trans = "log",
    # breaks = c(20, 50, 150),
    labels=scaleFUN) +
  coord_cartesian(ylim = c(-100, 50)) +
  scale_color_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")))) +
  scale_fill_manual(values = alpha(mod_cols, alpha = 0.5), labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = 6),
        plot.margin = margin(t = 16,  # Top margin
                             r = 123,  # Right margin
                             b = 8,  # Bottom margin
                             l = 11)  # Left margin
  ) +
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 50"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p2 <- ggplot(data = newdat %>% filter(scen_pa == 200), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s)))) +
  scale_y_continuous(
    # limits = ranges$`200`,
    name = expression(bold(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")"))),
    # trans = "log",
    # breaks = c(20, 50, 150),
    labels=scaleFUN) +
  coord_cartesian(ylim = c(-180, 90)) +
  scale_color_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")))) +
  scale_fill_manual(values = alpha(model_cols, alpha = 0.5), labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, axis.title.y = element_text(vjust = 6),
        plot.margin = margin(t = 8,  # Top margin
                             r = 0,  # Right margin
                             b = 16,  # Bottom margin
                             l = 11)  # Left margin
  ) + 
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 200"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p3 <- ggplot(data = newdat %>% filter(scen_pa == 800), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s)))) +
  scale_y_continuous(
    # limits = ranges$`800`,
    name = expression(bold(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")"))),
    # trans = "log",
    # breaks = c(20, 50, 150),
    labels=scaleFUN) +
  coord_cartesian(ylim = c(-26000, 14000)) +
  scale_color_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")))) +
  scale_fill_manual(values = alpha(model_cols, alpha = 0.5), labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none",
        plot.margin = margin(t = 0,  # Top margin
                             r = 123,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 800"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

library(gridExtra)
png(filename = paste0(home.wd, "/Figures/our_sim_results_boxplots_append.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
grid.arrange(p1,p2,p3)
dev.off()

# Original

# set the ranges for each n[PA]
tmp <- boxplot(KLdiv ~ scen_pa, data = plotdat)
ranges <- data.frame(apply(tmp$stats, 2, function(x){c(floor(x[1]), ceiling(x[5]))}))
colnames(ranges) <- c("50", "200", "800")

p1 <- ggplot(data = newdat %>% filter(scen_pa == 50), aes(x = b_range, y = KLdiv, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s)))) +
  scale_y_continuous(limits = ranges$`50`,
                     name = expression(bold(paste(D[KL],"(", mu, " || ", hat(mu), ")"))),
                     trans = "log",
                     labels=scaleFUN,
                     breaks = c(20, 50, 150)) +
  scale_color_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  scale_fill_manual(values = alpha(model_cols, alpha = 0.5), labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = 6),
        plot.margin = margin(t = 16,  # Top margin
                             r = 123,  # Right margin
                             b = 8,  # Bottom margin
                             l = 11)  # Left margin
  ) +
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 50")))

p2 <- ggplot(data = plotdat %>% filter(scen_pa == 200), aes(x = b_range, y = KLdiv, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s)))) +
  scale_y_continuous(limits = ranges$`200`,
                     name = expression(bold(paste(D[KL],"(", mu, " || ", hat(mu), ")"))),
                     trans = "log",
                     labels=scaleFUN,
                     breaks = c(50, 150, 300)) +
  scale_color_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  scale_fill_manual(values = alpha(model_cols, alpha = 0.5), labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, axis.title.y = element_text(vjust = 6),
        plot.margin = margin(t = 8,  # Top margin
                             r = 0,  # Right margin
                             b = 16,  # Bottom margin
                             l = 11)  # Left margin
  ) + 
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 200")))

p3 <- ggplot(data = plotdat %>% filter(scen_pa == 800), aes(x = b_range, y = KLdiv, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s)))) +
  scale_y_continuous(limits = ranges$`800`,
                     name = expression(bold(paste(D[KL],"(", mu, " || ", hat(mu), ")"))),
                     trans = "log",
                     labels=scaleFUN,
                     breaks = c(500, 2000, 5000, 20000)) +
  scale_color_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  scale_fill_manual(values = alpha(model_cols, alpha = 0.5), labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none",
        plot.margin = margin(t = 0,  # Top margin
                             r = 123,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 800")))

library(gridExtra)
png(filename = paste0(home.wd, "/Figures/our_sim_results_boxplots.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
grid.arrange(p1,p2,p3)
dev.off()

png(filename = paste0(home.wd, "/Figures/our_sim_results_boxplots.png"), width = 6.3 * plot.res, height = 3.5 * plot.res, res = plot.res)
ggplot(data = plotdat, aes(x = bias_range, y = KLdiv, col = Model)) +
  geom_boxplot(outlier.shape = NA) + facet_wrap(~ scen_pa, scales = "free_y", labeller = label_bquote(E*"["*n[PA]*"]"==.(as.numeric(as.character(scen_pa))))) +
  scale_x_continuous(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s)))) +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))
  scale_y_continuous(limits = c(0, 500), name = expression(bold(paste(D[KL],"(", mu, " || ", hat(mu), ")"))), breaks = c(seq(55, 80, by = 5), seq(140, 190, by = 10), seq(4e3, 10e3, by = 2e3), 15e3, 20e3), trans = "log", labels=scaleFUN) +
  scale_color_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  scale_fill_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10),
        legend.text=element_text(size=10), legend.text.align = 0, plot.title = element_text(hjust = 0.5))
dev.off()

# Old Plot of realised presence rate
dat$pres_rate <- dat$n_po / dat$n_pa_pres
dat$pres_rate_f <- cut(dat$pres_rate, breaks = c(0, 0.05, 0.25, 1, 4, 16, max(dat$pres_rate)))
levels(dat$pres_rate_f) <- c("< 0.05", "= (0.05, 0.25]", "= (0.25, 1]", "= (1, 4]", "= (4, 16]", "> 16")
model_cols <- c("darkgoldenrod1", "darkblue")
kl_breaks <- c(seq(0, 400, by = 30), seq(700, 1000, by = 50), seq(2000, 8000, by = 200))
png(filename = paste0(home.wd, "/Figures/our_sim_results_pres_ratio.png"), width = 6.3 * plot.res, height = 4.5 * plot.res, res = plot.res)
ggplot(data = dat %>% filter(Model != "PO Only"), aes(x = bias_range, y = KLdiv, col = Model)) +
  geom_smooth(aes(fill = Model)) + facet_wrap(~ pres_rate_f, scales = "free_y", labeller = label_bquote(frac(n[PO], Sigma~p[i])~.(as.character(pres_rate_f)))) +
  scale_x_continuous(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s)))) +
  scale_y_continuous(name = expression(bold(paste(D[KL],"(", mu, " || ", hat(mu), ")"))), breaks = kl_breaks, trans = "log", labels=scaleFUN) +
  scale_color_manual(values = model_cols) + scale_fill_manual(values = model_cols) +
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
              axis.text.y = element_text(size=10), axis.text.x = element_text(size=10),
              legend.text=element_text(size=10), legend.text.align = 0, plot.title = element_text(hjust = 0.5))
dev.off()

# the ratios considered in simulations
sort(unique(apply(expand.grid(c(50, 200, 800), c(50, 200, 800)), 1, function(x){x[1]/x[2]})))

################################################################################
# Append Full Results Table ####################################################
################################################################################

tab <- dat %>% filter(Model != "PO Only") %>% group_by(scen, env_range, latent_range, Model) %>%
  summarise(KL = median(KLdiv), MAE = median(MAE), Comp_Time = median(timing))
# set up a unique id to link to fit issues table
id1 <- paste(tab$scen, tab$env_range, tab$latent_range, tab$Model, sep = ":")
id2 <- paste(fit.issues$scen, fit.issues$env_range, fit.issues$latent_range, fit.issues$Model, sep = ":")
# add in the fit issues
tab <- cbind(tab, fit.issues[match(id1, id2), c("N_SIMS", "INF_KL", "POOR_CONV")])

# separate the expected sample size scenarios
tab$n_pa <- sapply(strsplit(tab$scen, "PA", fixed = T), function(x){x[2]})
tab$n_po <- sapply(strsplit(sapply(strsplit(tab$scen, "_", fixed = T), function(x){x[1]}), "PO", fixed = T), function(x){x[2]})

# re-order and align the columns for presenting
tab <- tab[order(tab$n_pa, tab$n_po, decreasing = T), c("n_pa", "n_po", "env_range", "latent_range", "Model", "KL", "Comp_Time", "INF_KL", "POOR_CONV")]

library(xtable)
print(xtable(tab, label = "append:tab:sim", digits = 2,
             caption = c("Full results from our simulations")),
      include.rownames=F, caption.placement = "top")

################################################################################
# Append plot: Difference in shared effect ranges ##############################
################################################################################

# adjust the plotting data for nice labels within facet_grid() for both pa_scen and env_range
plotdat2 <- dat[dat$Model != "PO Only", ]
plotdat2$pa_scen <- factor(plotdat2$scen_pa,
                          levels = c("50", "200", "800"),
                          labels = c("E*'['*n[PA]*']'==50", "E*'['*n[PA]*']'==200", "E*'['*n[PA]*']'==800")
)
plotdat2$env_range <- factor(plotdat2$env_range,
                             levels = c(15, 25),
                             labels = c("rho[X]==15~~(Env.)~~~~rho[xi]==25~~(Latent)",
                                        "rho[X]==25~~(Env.)~~~~rho[xi]==15~~(Latent)"
                             )
)
# re-factor the model to incorporate PO scenario into IDM
plotdat2$Model[plotdat2$Model == "IDM"] <- paste(plotdat2$Model[plotdat2$Model == "IDM"], plotdat2$scen_po[plotdat2$Model == "IDM"], sep = " - ")
plotdat2$Model <- factor(plotdat2$Model,
                        levels = c("IDM - 50", "IDM - 200", "IDM - 800", "PA Only"),
                        labels =  c("IDM - N PO = 50", "IDM - N PO = 200", "IDM - N PO = 800", "PA only")
)

# set up the vertical lines for plotting to indicate share effect ranges
verticals <- rbind( 
  data.frame(intercepts=c(25,15), type=c("Latent", "Env."), scen = "N PO = 50\nN PA = 50", env_range = "Fine scale Environ. - Coarse scale Latent"),
  data.frame(intercepts=c(15,25), type=c("Latent", "Env."), scen = "N PO = 50\nN PA = 50", env_range = "Coarse scale Environ. - Fine scale Latent"),
  data.frame(intercepts=c(25,15), type=c("Latent", "Env."), scen = "N PO = 200\nN PA = 50", env_range = "Fine scale Environ. - Coarse scale Latent"),
  data.frame(intercepts=c(15,25), type=c("Latent", "Env."), scen = "N PO = 200\nN PA = 50", env_range = "Coarse scale Environ. - Fine scale Latent"),
  data.frame(intercepts=c(25,15), type=c("Latent", "Env."), scen = "N PO = 800\nN PA = 50", env_range = "Fine scale Environ. - Coarse scale Latent"),
  data.frame(intercepts=c(15,25), type=c("Latent", "Env."), scen = "N PO = 800\nN PA = 50", env_range = "Coarse scale Environ. - Fine scale Latent"),
  data.frame(intercepts=c(25,15), type=c("Latent", "Env."), scen = "N PO = 50\nN PA = 200", env_range = "Fine scale Environ. - Coarse scale Latent"),
  data.frame(intercepts=c(15,25), type=c("Latent", "Env."), scen = "N PO = 50\nN PA = 200", env_range = "Coarse scale Environ. - Fine scale Latent"),
  data.frame(intercepts=c(25,15), type=c("Latent", "Env."), scen = "N PO = 200\nN PA = 200", env_range = "Fine scale Environ. - Coarse scale Latent"),
  data.frame(intercepts=c(15,25), type=c("Latent", "Env."), scen = "N PO = 200\nN PA = 200", env_range = "Coarse scale Environ. - Fine scale Latent"),
  data.frame(intercepts=c(25,15), type=c("Latent", "Env."), scen = "N PO = 800\nN PA = 200", env_range = "Fine scale Environ. - Coarse scale Latent"),
  data.frame(intercepts=c(15,25), type=c("Latent", "Env."), scen = "N PO = 800\nN PA = 200", env_range = "Coarse scale Environ. - Fine scale Latent"),
  data.frame(intercepts=c(25,15), type=c("Latent", "Env."), scen = "N PO = 50\nN PA = 800", env_range = "Fine scale Environ. - Coarse scale Latent"),
  data.frame(intercepts=c(15,25), type=c("Latent", "Env."), scen = "N PO = 50\nN PA = 800", env_range = "Coarse scale Environ. - Fine scale Latent"),
  data.frame(intercepts=c(25,15), type=c("Latent", "Env."), scen = "N PO = 200\nN PA = 800", env_range = "Fine scale Environ. - Coarse scale Latent"),
  data.frame(intercepts=c(15,25), type=c("Latent", "Env."), scen = "N PO = 200\nN PA = 800", env_range = "Coarse scale Environ. - Fine scale Latent"),
  data.frame(intercepts=c(25,15), type=c("Latent", "Env."), scen = "N PO = 800\nN PA = 800", env_range = "Fine scale Environ. - Coarse scale Latent"),
  data.frame(intercepts=c(15,25), type=c("Latent", "Env."), scen = "N PO = 800\nN PA = 800", env_range = "Coarse scale Environ. - Fine scale Latent")
)
verticals$pa_scen <- factor(unlist(lapply(strsplit(as.character(verticals$scen), "\n", fixed = T), function(x){x[2]})),
                           levels = c("N PA = 50", "N PA = 200", "N PA = 800"),
                           labels = c("E*'['*n[PA]*']'==50", "E*'['*n[PA]*']'==200", "E*'['*n[PA]*']'==800")
)
verticals$env_range <- factor(verticals$env_range,
                             levels = c("Fine scale Environ. - Coarse scale Latent",
                                        "Coarse scale Environ. - Fine scale Latent"
                             ),
                             labels = c("rho[X]==15~~(Env.)~~~~rho[xi]==25~~(Latent)",
                                        "rho[X]==25~~(Env.)~~~~rho[xi]==15~~(Latent)"
                             )
)
verticals$type <- factor(verticals$type,
                                   levels = c("Env.", "Latent"),
                                   labels = c("Env.~(X)", "Latent~(xi)")
)
model_cols <- c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4", "darkblue")

png(filename = paste0(home.wd, "/Figures/our_sim_results_fixed_effect_range.png"), width = 6.3 * plot.res, height = 5 * plot.res, res = plot.res)
ggplot(data = plotdat2, aes(x = bias_range, y = KLdiv, col = Model)) +
  geom_smooth(aes(fill = Model)) + facet_grid(vars(env_range), vars(pa_scen), labeller = label_parsed) +
  scale_x_continuous(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s)))) +
  scale_y_continuous(name = expression(paste(D[KL],"(", mu, " || ", hat(mu), ")")), breaks = c(100, 1000, 10000, 100000), trans = "log", labels=scaleFUN) +
  scale_color_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  scale_fill_manual(values = model_cols, labels = c(expression(paste("IDM: ", E*"["*n[PO]*"]", " = 50")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 200")), expression(paste("IDM: ", E*"["*n[PO]*"]", " = 800")), expression("PA only"))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10),
        legend.text=element_text(size=10), legend.text.align = 0, plot.title = element_text(hjust = 0.5)) +
  geom_vline(
    data=verticals,
    mapping=aes(xintercept=intercepts, linetype=type),
    size=0.8, color="black", #linetype = c("dashed", "dotted"),
    key_glyph="path"   # this makes the legend key horizontal lines, not vertical
  ) + scale_linetype_manual(name = "Shared\nEffect\nRange", values = c("dashed", "dotted"), labels = c("Env. (X)", expression(paste("Latent (", xi, ")")))) #+ ggtitle("PA Data Presences")
dev.off()