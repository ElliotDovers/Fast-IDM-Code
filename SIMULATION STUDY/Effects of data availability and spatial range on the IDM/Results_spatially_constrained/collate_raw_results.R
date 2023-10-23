## Analyse Scenario ##

home.wd <- getwd()

# # initialise the result storage
# dat <- NULL
# 
# # change into appropriate result folder
# setwd(paste(home.wd, "raw", sep = "/"))
# # get a list of all the individual result files
# res.list <- list.files()[grepl("res_", list.files(), fixed = T)]
# # inner loop through individual sim-by-model files
# for (job in res.list) {
#   # get the job number
#   job.no <- as.numeric(unlist(lapply(strsplit(unlist(lapply(strsplit(job, ".", fixed = T), function(x){x[1]})), "_", fixed = T), function(x){x[length(x)]})))
#   # load in the individual simulation results
#   load(job)
#   # add in the job number
#   res_tab$job <- job.no
#   # add in the scenario
#   res_tab$scen <- paste0("PO", res_tab$scen_po, "_PA", res_tab$scen_pa)
#   # add to the data
#   dat <- rbind(dat, res_tab)
#   # ranges <- rbind(ranges, attr(res_tab, "correlation ranges"))
#   rm(res_tab)
# }
# 
# setwd(home.wd)
# 
# save(list = "dat", file = "collated raw results.RDATA")
load("collated raw results.RDATA")

library(dplyr)
library(ggplot2)

dat$scenario <- factor(paste0(dat$scen, "_", dat$env_range))

# check out the number of flagged convergences - this are fit failures for scampr
dat %>% group_by(scen, Model, env_range, latent_range) %>%
  summarise(FLAGGED = sum(fit_flag == 1), INF_KL_sc = sum(is.infinite(KLdiv_sc)), INF_KL_po = sum(is.infinite(KLdiv_po)), INF_LL = sum(is.infinite(LL)), POOR_CONV = sum(LL > 1e5 | is.nan(LL)))

# create full results table for Appendix
fit.issues <- dat %>% group_by(scen, env_range, latent_range, Model) %>%
  summarise(N_SIMS = length((sim)), INF_KL_sc = sum(is.infinite(KLdiv_sc) & !(LL > 1e5 | is.nan(LL))), POOR_CONV = sum(LL > 1e5 | is.nan(LL)))

# check that the infinite KLs are not concentrated so that there are not enough simulations
sim.check <- dat %>% filter(is.finite(KLdiv_sc)) %>% group_by(Model, scen, env_range, bias_range) %>%
  summarise(unique_sims = length(unique(sim)), sims = length(sim), avg_k = mean(k), med_k = median(k), avg_k_bias = mean(k_bias), med_k_bias = median(k_bias))
min(sim.check$sims) # 95 is the lowest number returned in any one combination of env, lat and bias ranges
sim.check <- dat %>% filter(is.finite(KLdiv_po)) %>% group_by(Model, scen, env_range, bias_range) %>%
  summarise(unique_sims = length(unique(sim)), sims = length(sim), avg_k = mean(k), med_k = median(k), avg_k_bias = mean(k_bias), med_k_bias = median(k_bias))
min(sim.check$sims) # 95 is the lowest number returned in any one combination of env, lat and bias ranges

# check the number of basis functions selected for the problem fits
dat %>% filter(is.finite(KLdiv_sc) | is.na(KLdiv_sc)) %>% group_by(Model, scen, env_range, bias_range) %>%
  summarise(unique_sims = length(unique(sim)), sims = length(sim), avg_k = mean(k), med_k = median(k), avg_k_bias = mean(k_bias), med_k_bias = median(k_bias))
# the majority of problem KL are due to poor convergences
table(dat[(is.infinite(dat$KLdiv_sc) | is.na(dat$KLdiv_sc)) & dat$Model %in% c("PA Only SC", "IDM SC"), "LL"] > 1e5)
table(dat[(is.infinite(dat$KLdiv_po) | is.na(dat$KLdiv_po)) & dat$Model %in% c("PO Only", "IDM SC"), "LL"] > 1e5)

dat$id <- paste(dat$scen, dat$env_range, dat$bias_range, dat$sim, sep = ":")
# remove entire sims with infinite or NaN KLdiv (SO TO COMPARE MODELS FAIRLY)
tmp <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(check = sum(is.infinite(KLdiv_full) | is.na(KLdiv_full)))
rm.ids <- paste(tmp$scen[tmp$check != 0], tmp$env_range[tmp$check != 0], tmp$bias_range[tmp$check != 0], tmp$sim[tmp$check != 0], sep = ":")
tmp <- dat %>% filter(Model %in% c("PA Only SC", "IDM SC")) %>% group_by(scen, env_range, bias_range, sim) %>% summarise(check = sum(is.infinite(KLdiv_sc) | is.na(KLdiv_sc)))
rm.ids_sc <- paste(tmp$scen[tmp$check != 0], tmp$env_range[tmp$check != 0], tmp$bias_range[tmp$check != 0], tmp$sim[tmp$check != 0], sep = ":")
tmp <- dat %>% filter(Model %in% c("PO Only", "IDM SC")) %>% group_by(scen, env_range, bias_range, sim) %>% summarise(check = sum(is.infinite(KLdiv_po) | is.na(KLdiv_po)))
rm.ids_po <- paste(tmp$scen[tmp$check != 0], tmp$env_range[tmp$check != 0], tmp$bias_range[tmp$check != 0], tmp$sim[tmp$check != 0], sep = ":")
identical(rm.ids_sc, rm.ids_po) # all the failed sims are the same for either comparison we wish to make
all(rm.ids_po %in% rm.ids) # and all of these are also failures for the full KL metrics
dat <- dat[!dat$id %in% rm.ids, ]

# take out only the models we wish to compare
dat <- dat[dat$Model %in% c("PA Only SC", "PO Only", "IDM SC"), ]

# there are no crazy large KL div from the spatially constrained models
dat[dat$KLdiv_full > 1e6 & dat$Model != "PO Only", ]
dat[dat$KLdiv_sc > 1e6 & dat$Model != "PO Only", ]
dat[dat$KLdiv_po > 1e6 & dat$Model != "PA Only SC", ]

# factor the environment variable range with nice labels
dat$fixed_ranges <- factor(dat$env_range,
                      levels = c("15",
                                 "25"
                      ),
                      labels = c("Fine scale Environ. - Coarse scale Latent",
                                 "Coarse scale Environ. - Fine scale Latent"
                      )
)

# adjust our metric to be relative to the PA Only model for mu(s)
dat_pa <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(delta_kl = KLdiv_full[Model == "IDM SC"] - KLdiv_full[Model == "PA Only SC"])
dat_pa$scen_pa <- sapply(strsplit(dat_pa$scen, "_PA", fixed = T), function(x){x[2]})
dat_pa$Model <- paste0("IDM - N PO = ", sapply(strsplit(sapply(strsplit(dat_pa$scen, "PO", fixed = T), function(x){x[2]}), "_", fixed = T), function(x){x[1]}))
dat_pa$b_range <- factor(dat_pa$bias_range)
dat_pa$Model <- factor(dat_pa$Model, levels = c("IDM - N PO = 50", "IDM - N PO = 200", "IDM - N PO = 800"))
dat_pa$scen_pa <- factor(dat_pa$scen_pa, levels = c(50, 200, 800))

# adjust our metric to be relative to the PA Only model for lambda(s)
dat_po <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(delta_kl = KLdiv_po[Model == "IDM SC"] - KLdiv_po[Model == "PO Only"])
dat_po$scen_po <- sapply(strsplit(sapply(strsplit(dat_po$scen, "PO", fixed = T), function(x){x[2]}), "_", fixed = T), function(x){x[1]})
dat_po$Model <- paste0("IDM - N PA = ", sapply(strsplit(dat_po$scen, "_PA", fixed = T), function(x){x[2]}))
dat_po$b_range <- factor(dat_po$bias_range)
dat_po$Model <- factor(dat_po$Model, levels = c("IDM - N PA = 50", "IDM - N PA = 200", "IDM - N PA = 800"))
dat_po$scen_po <- factor(dat_po$scen_po, levels = c(50, 200, 800))

# set colours
mod_cols_pa <- c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4")
mod_cols_po <- c("darkorchid2", "darkorchid3", "darkorchid4")
# set figure size
plot.res <- 500
# set a scale function to round to KL div. axis values to zero decimal places
scaleFUN <- function(x) sprintf("%.0f", x)

# set the ranges for each n[PA] and n[PO]
tmp <- boxplot(delta_kl ~ scen_pa, data = dat_pa)
pa_ranges <- data.frame(apply(tmp$stats, 2, function(x){c(floor(x[1]) + 0.05 * x[1], ceiling(x[5]) + 0.05 * x[5])}))
colnames(pa_ranges) <- c("pa50", "pa200", "pa800")
tmp <- boxplot(delta_kl ~ scen_po, data = dat_po)
po_ranges <- data.frame(apply(tmp$stats, 2, function(x){c(floor(x[1]) + 0.05 * x[1], ceiling(x[5]) + 0.05 * x[5])}))
colnames(po_ranges) <- c("po50", "po200", "po800")


Ap1 <- ggplot(data = dat_pa %>% filter(scen_pa == 50 & env_range == 15), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s)))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = pa_ranges$pa50) +
  scale_color_manual(name = "Data\nIntegrated", values = mod_cols_pa, labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_fill_manual(name = "Data\nIntegrated", values = alpha(mod_cols_pa, alpha = 0.5), labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = 0),
        plot.margin = margin(t = 16,  # Top margin
                             r = 88,  # Right margin
                             b = 8,  # Bottom margin
                             l = 0)  # Left margin
  ) +
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 50"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

Ap2 <- ggplot(data = dat_pa %>% filter(scen_pa == 200 & env_range == 15), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s)))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = pa_ranges$pa200) +
  scale_color_manual(name = "Data\nIntegrated", values = mod_cols_pa, labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_fill_manual(name = "Data\nIntegrated", values = alpha(mod_cols_pa, alpha = 0.5), labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=8), legend.text.align = 0, axis.title.y = element_text(vjust = 0),
        plot.margin = margin(t = 8,  # Top margin
                             r = 0,  # Right margin
                             b = 16,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 200"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

Ap3 <- ggplot(data = dat_pa %>% filter(scen_pa == 800 & env_range == 15), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s))) , breaks = seq(11, 29, by = 4)) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")), breaks = c(-20000, -10000, 0, 10000), labels = c("-20K", "-10K", "0", "10K")) +
  coord_cartesian(ylim = pa_ranges$pa800) +
  scale_color_manual(name = "Data\nIntegrated", values = mod_cols_pa, labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_fill_manual(name = "Data\nIntegrated", values = alpha(mod_cols_pa, alpha = 0.5), labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none",
        plot.margin = margin(t = 0,  # Top margin
                             r = 88,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 800"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

Ap4 <- ggplot(data = dat_po %>% filter(scen_po == 50 & env_range == 15), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s)))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = po_ranges$po50) +
  scale_color_manual(name = "Data\nIntegrated", values = mod_cols_po, labels = c(expression(paste(E*"["*n[PA]*"]", " = 50")), expression(paste(E*"["*n[PA]*"]", " = 200")), expression(paste(E*"["*n[PA]*"]", " = 800")))) +
  scale_fill_manual(name = "Data\nIntegrated", values = alpha(mod_cols_po, alpha = 0.5), labels = c(expression(paste(E*"["*n[PA]*"]", " = 50")), expression(paste(E*"["*n[PA]*"]", " = 200")), expression(paste(E*"["*n[PA]*"]", " = 800")))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = -1),
        plot.margin = margin(t = 16,  # Top margin
                             r = 88,  # Right margin
                             b = 8,  # Bottom margin
                             l = 0)  # Left margin
  ) +
  ggtitle(label = expression(paste(E*"["*n[PO]*"]", " = 50"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

Ap5 <- ggplot(data = dat_po %>% filter(scen_po == 200 & env_range == 15), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s)))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = po_ranges$po200) +
  scale_color_manual(name = "Data\nIntegrated", values = mod_cols_po, labels = c(expression(paste(E*"["*n[PA]*"]", " = 50")), expression(paste(E*"["*n[PA]*"]", " = 200")), expression(paste(E*"["*n[PA]*"]", " = 800")))) +
  scale_fill_manual(name = "Data\nIntegrated", values = alpha(mod_cols_po, alpha = 0.5), labels = c(expression(paste(E*"["*n[PA]*"]", " = 50")), expression(paste(E*"["*n[PA]*"]", " = 200")), expression(paste(E*"["*n[PA]*"]", " = 800")))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=8), legend.text.align = 0, axis.title.y = element_text(vjust = -1),
        plot.margin = margin(t = 8,  # Top margin
                             r = 0,  # Right margin
                             b = 16,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  ggtitle(label = expression(paste(E*"["*n[PO]*"]", " = 200"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

Ap6 <- ggplot(data = dat_po %>% filter(scen_po == 800 & env_range == 15), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s))) , breaks = seq(11, 29, by = 4)) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = po_ranges$po800) +
  scale_color_manual(name = "Data\nIntegrated", values = mod_cols_po, labels = c(expression(paste(E*"["*n[PA]*"]", " = 50")), expression(paste(E*"["*n[PA]*"]", " = 200")), expression(paste(E*"["*n[PA]*"]", " = 800")))) +
  scale_fill_manual(name = "Data\nIntegrated", values = alpha(mod_cols_po, alpha = 0.5), labels = c(expression(paste(E*"["*n[PA]*"]", " = 50")), expression(paste(E*"["*n[PA]*"]", " = 200")), expression(paste(E*"["*n[PA]*"]", " = 800")))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=10), axis.title.y = element_text(vjust = -1),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none",
        plot.margin = margin(t = 0,  # Top margin
                             r = 88,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  ggtitle(label = expression(paste(E*"["*n[PO]*"]", " = 800"))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

library(gridExtra)
png(filename = paste0(home.wd, "/Figures/our_sim_results_boxplots_append_short_env_spatially_constrained.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
grid.arrange(Ap1,Ap4,Ap2,Ap5,Ap3,Ap6, nrow = 3, ncol = 2)
dev.off()

# adjust our metric to be relative to the PA Only model
dat_pa_idm <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(delta_kl = KLdiv_full[Model == "IDM SC"] - KLdiv_full[Model == "PA Only SC"])
# dat_pa_idm <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(delta_kl = KLdiv_sc[Model == "IDM SC"] - KLdiv_sc[Model == "PA Only SC"])
# dat_pa_idm$scen_pa <- sapply(strsplit(dat_pa_idm$scen, "_PA", fixed = T), function(x){x[2]})
# dat_pa_idm$Model <- paste0("IDM - N PO = ", sapply(strsplit(sapply(strsplit(dat_pa_idm$scen, "PO", fixed = T), function(x){x[2]}), "_", fixed = T), function(x){x[1]}))
# dat_pa_idm$b_range <- factor(dat_pa_idm$bias_range)
# dat_pa_idm$Model <- factor(dat_pa_idm$Model, levels = c("IDM - N PO = 50", "IDM - N PO = 200", "IDM - N PO = 800"))
dat_pa_idm <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(delta_kl = KLdiv_full[Model == "IDM SC"] - KLdiv_full[Model == "PA Only SC"])
dat_pa_idm$scen_pa <- sapply(strsplit(dat_pa_idm$scen, "_PA", fixed = T), function(x){x[2]})
dat_pa_idm$Model <- paste0("IDM - N PO = ", sapply(strsplit(sapply(strsplit(dat_pa_idm$scen, "PO", fixed = T), function(x){x[2]}), "_", fixed = T), function(x){x[1]}))
dat_pa_idm$b_range <- factor(dat_pa_idm$bias_range)
dat_pa_idm$Model <- factor(dat_pa_idm$Model, levels = c("IDM - N PO = 50", "IDM - N PO = 200", "IDM - N PO = 800"))

# set the plot resolution
plot.res <- 500

# set the ranges for each n[PA]
tmp <- boxplot(delta_kl ~ scen_pa, data = dat_pa_idm, plot = F)
ranges <- data.frame(apply(tmp$stats, 2, function(x){c(floor(x[1]), ceiling(x[5]))}))
colnames(ranges) <- c("50", "200", "800")
dat_pa_idm$scen_pa <- factor(dat_pa_idm$scen_pa, levels = c(50, 200, 800))

mod_cols <- c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4")

p1 <- ggplot(data = dat_pa_idm %>% filter(scen_pa == 50 & b_range == "19"), aes(x = Model, y = delta_kl, color = Model, fill = Model)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected Amount of Presence-Only Data Integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_y_continuous(
    # limits = ranges$`200`,
    name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")),
    # trans = "log",
    # breaks = c(20, 50, 150),
    labels=scaleFUN) +
  coord_cartesian(ylim = c(-480, 280)) +
  scale_color_manual(values = mod_cols) +
  scale_fill_manual(values = alpha(mod_cols, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = 6),
        plot.margin = margin(t = 20,  # Top margin
                             r = 15,  # Right margin
                             b = 12,  # Bottom margin
                             l = 12)  # Left margin
  ) +
  # ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 50"))) +
  ggtitle(label = "Species expected at 5% of survey sites") +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p2 <- ggplot(data = dat_pa_idm %>% filter(scen_pa == 200 & b_range == "19"), aes(x = Model, y = delta_kl, color = Model, fill = Model)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected Amount of Presence-Only Data Integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_y_continuous(
    # limits = ranges$`200`,
    name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")),
    # trans = "log",
    # breaks = c(20, 50, 150),
    labels=scaleFUN) +
  coord_cartesian(ylim = c(-2000, 1000)) +
  scale_color_manual(values = mod_cols) +
  scale_fill_manual(values = alpha(mod_cols, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = 3),
        plot.margin = margin(t = 12,  # Top margin
                             r = 15,  # Right margin
                             b = 20,  # Bottom margin
                             l = 5)  # Left margin
  ) +
  # ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 200"))) +
  ggtitle(label = "Species expected at 20% of survey sites") +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p3 <- ggplot(data = dat_pa_idm %>% filter(scen_pa == 800 & b_range == "19"), aes(x = Model, y = delta_kl, col = Model, fill = Model)) +
  # geom_boxplot(outlier.shape = NA) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected Amount of Presence-Only Data Integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_y_continuous(
    # limits = ranges$`800`,
    name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")),
    # trans = "log",
    # breaks = c(20, 50, 150),
    labels=scaleFUN) +
  coord_cartesian(ylim = c(-40000, 25000)) +
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
  # ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 800"))) +
  ggtitle(label = "Species expected at 80% of survey sites") +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

library(gridExtra)
png(filename = paste0(getwd(), "/Figures/sim_results_spatially_constrained.png"), width = 4.2 * plot.res, height = 6.5 * plot.res, res = plot.res)
grid.arrange(p1,p2,p3)
dev.off()

# check there are no crazy large values (indicating poor convergence)
dat[dat$KLdiv_full > 1e6 & !dat$Model %in% c("PO Only"), ]

dat_pa_idm <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(sc = KLdiv_full[Model == "IDM SC"] - KLdiv_full[Model == "PA Only SC"], full = KLdiv_full[Model == "IDM FULL"] - KLdiv_full[Model == "PA Only FULL"])
plotdat_pa <- data.table::melt(data = dat_pa_idm, measure.vars = c("sc", "full"))
plotdat_pa$brange <- factor(plotdat_pa$bias_range)
plotdat_pa$scen_po <- factor(sapply(strsplit(sapply(strsplit(plotdat_pa$scen, "_PA", fixed = T), function(x){x[1]}), "PO", fixed = T), function(x){x[2]}),
                             levels = c("50", "200", "800"),
                             labels = c("E*'['*n[PO]*']'==50", "E*'['*n[PO]*']'==200", "E*'['*n[PO]*']'==800")
)
plotdat_pa$scen_pa <- factor(sapply(strsplit(plotdat_pa$scen, "_PA", fixed = T), function(x){x[2]}),
                             levels = c("50", "200", "800"),
                             labels = c("E*'['*n[PA]*']'==50", "E*'['*n[PA]*']'==200", "E*'['*n[PA]*']'==800")
)
plot.res <- 500
mod_cols <- c("darkorange", "royalblue")
png(filename = paste0(getwd(), "/Figures/sim_results_spatially_constrained_lambda.png"), width = 6.3 * plot.res, height = 5.5 * plot.res, res = plot.res)
ggplot(data = plotdat_pa %>% filter(env_range == 15# & bias_range == 19
), aes(x = brange, y = value, color = variable, fill = variable)) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  geom_boxplot(outlier.shape = NA) +
  # facet_wrap(~ scen, scales = "free")#, var(brange))
  facet_grid(vars(scen_po), vars(scen_pa), labeller = label_parsed) +
  coord_cartesian(ylim = c(-1000, 500)) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s))), breaks = seq(11, 29, by = 4)) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")"))) +
  scale_color_manual(name = "Spatially\nConstrained\nPA Data", values = mod_cols, labels = c("Yes", "No")) +
  scale_fill_manual(name = "Spatially\nConstrained\nPA Data", values = alpha(mod_cols, alpha = 0.5), labels = c("Yes", "No")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10)
  ) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")
dev.off()


## Additionally look at improvement in the PO intensity estimation #############

# check there are no crazy large values (indicating poor convergence)
dat[dat$KLdiv_po > 1e6 & !dat$Model %in% c("PA Only SC", "PA Only FULL"), ]

dat_po_idm <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(sc = KLdiv_po[Model == "IDM SC"] - KLdiv_po[Model == "PO Only"], full = KLdiv_po[Model == "IDM FULL"] - KLdiv_po[Model == "PO Only"])
plotdat_po <- data.table::melt(data = dat_po_idm, measure.vars = c("sc", "full"))
plotdat_po$brange <- factor(plotdat_po$bias_range)
plotdat_po$scen_po <- factor(sapply(strsplit(sapply(strsplit(plotdat_po$scen, "_PA", fixed = T), function(x){x[1]}), "PO", fixed = T), function(x){x[2]}),
                             levels = c("50", "200", "800"),
                             labels = c("E*'['*n[PO]*']'==50", "E*'['*n[PO]*']'==200", "E*'['*n[PO]*']'==800")
)
plotdat_po$scen_pa <- factor(sapply(strsplit(plotdat_po$scen, "_PA", fixed = T), function(x){x[2]}),
                             levels = c("50", "200", "800"),
                             labels = c("E*'['*n[PA]*']'==50", "E*'['*n[PA]*']'==200", "E*'['*n[PA]*']'==800")
)
plot.res <- 500
mod_cols <- c("darkorange", "royalblue")
png(filename = paste0(getwd(), "/Figures/sim_results_spatially_constrained_lambda.png"), width = 6.3 * plot.res, height = 5.5 * plot.res, res = plot.res)
ggplot(data = plotdat_po %>% filter(env_range == 15# & bias_range == 19
                                    ), aes(x = brange, y = value, color = variable, fill = variable)) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  geom_boxplot(outlier.shape = NA) +
  # facet_wrap(~ scen, scales = "free")#, var(brange))
  facet_grid(vars(scen_po), vars(scen_pa), labeller = label_parsed) +
  coord_cartesian(ylim = c(-10, 5)) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", B(s))), breaks = seq(11, 29, by = 4)) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")"))) +
  scale_color_manual(name = "Spatially\nConstrained\nPA Data", values = mod_cols, labels = c("Yes", "No")) +
  scale_fill_manual(name = "Spatially\nConstrained\nPA Data", values = alpha(mod_cols, alpha = 0.5), labels = c("Yes", "No")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10)
  ) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")
dev.off()
# # OLD PLOT
# dat_po_idm <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(delta_kl = KLdiv_po[Model == "IDM SC"] - KLdiv_po[Model == "PO Only"])
# dat_po_idm$scen_po <- sapply(strsplit(sapply(strsplit(dat_po_idm$scen, "_PA", fixed = T), function(x){x[1]}), "PO", fixed = T), function(x){x[2]})
# dat_po_idm$Model <- paste0("IDM - N PA = ", sapply(strsplit(dat_po_idm$scen, "_PA", fixed = T), function(x){x[2]}))
# dat_po_idm$b_range <- factor(dat_po_idm$bias_range)
# dat_po_idm$Model <- factor(dat_po_idm$Model, levels = c("IDM - N PA = 50", "IDM - N PA = 200", "IDM - N PA = 800"))
# 
# # set the plot resolution
# plot.res <- 500
# 
# # set the ranges for each n[PA]
# tmp <- boxplot(delta_kl ~ scen_pa, data = dat_po_idm, plot = F)
# ranges <- data.frame(apply(tmp$stats, 2, function(x){c(floor(x[1]), ceiling(x[5]))}))
# colnames(ranges) <- c("50", "200", "800")
# dat_pa_idm$scen_pa <- factor(dat_pa_idm$scen_pa, levels = c(50, 200, 800))
# 
# mod_cols <- c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4")
# 
# p1 <- ggplot(data = dat_po_idm %>% filter(scen_po == 50 & b_range == "19"), aes(x = Model, y = delta_kl, color = Model, fill = Model)) +
#   # geom_boxplot(outlier.shape = NA) +
#   geom_boxplot(outlier.shape = 1, outlier.size = 1) +
#   scale_x_discrete(name = "Expected Amount of Presence-Only Data Integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
#   scale_y_continuous(
#     # limits = ranges$`200`,
#     name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")")),
#     # trans = "log",
#     # breaks = c(20, 50, 150),
#     labels=scaleFUN) +
#   coord_cartesian(ylim = c(-15, 10)) +
#   scale_color_manual(values = mod_cols) +
#   scale_fill_manual(values = alpha(mod_cols, alpha = 0.5)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
#         axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
#         legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", #axis.title.y = element_text(vjust = 6),
#         plot.margin = margin(t = 20,  # Top margin
#                              r = 15,  # Right margin
#                              b = 12,  # Bottom margin
#                              l = 0)  # Left margin
#   ) +
#   # ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 50"))) +
#   ggtitle(label = expression(paste(E*"["*n[PO]*"]", " = 50 expected presence-only data"))) +
#   geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")
# 
# p2 <- ggplot(data = dat_po_idm %>% filter(scen_po == 200 & b_range == "19"), aes(x = Model, y = delta_kl, color = Model, fill = Model)) +
#   # geom_boxplot(outlier.shape = NA) +
#   geom_boxplot(outlier.shape = 1, outlier.size = 1) +
#   scale_x_discrete(name = "Expected Amount of Presence-Only Data Integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
#   scale_y_continuous(
#     # limits = ranges$`200`,
#     name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")")),
#     # trans = "log",
#     # breaks = c(20, 50, 150),
#     labels=scaleFUN) +
#   coord_cartesian(ylim = c(-15, 10)) +
#   scale_color_manual(values = mod_cols) +
#   scale_fill_manual(values = alpha(mod_cols, alpha = 0.5)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
#         axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
#         legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", #axis.title.y = element_text(vjust = 6),
#         plot.margin = margin(t = 12,  # Top margin
#                              r = 15,  # Right margin
#                              b = 20,  # Bottom margin
#                              l = 0)  # Left margin
#   ) +
#   # ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 200"))) +
#   ggtitle(label = expression(paste(E*"["*n[PO]*"]", " = 200 expected presence-only data"))) +
#   geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")
# 
# p3 <- ggplot(data = dat_po_idm %>% filter(scen_po == 800 & b_range == "19"), aes(x = Model, y = delta_kl, col = Model, fill = Model)) +
#   # geom_boxplot(outlier.shape = NA) +
#   geom_boxplot(outlier.shape = 1, outlier.size = 1) +
#   scale_x_discrete(name = "Expected presences in the survey data", labels = c(expression(paste(E*"["*n[PA]*"]", " = 50")), expression(paste(E*"["*n[PA]*"]", " = 200")), expression(paste(E*"["*n[PA]*"]", " = 800")))) +
#   scale_y_continuous(
#     # limits = ranges$`800`,
#     name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")")),
#     # trans = "log",
#     # breaks = c(20, 50, 150),
#     labels=scaleFUN) +
#   coord_cartesian(ylim = c(-15, 10)) +
#   scale_color_manual(values = mod_cols) +
#   scale_fill_manual(values = alpha(mod_cols, alpha = 0.5)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
#         axis.text.y = element_text(size=10), axis.text.x = element_text(size=10), axis.title.x = element_text(margin = margin(t = 15)),
#         legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none",
#         plot.margin = margin(t = 0,  # Top margin
#                              r = 15,  # Right margin
#                              b = 0,  # Bottom margin
#                              l = 0)  # Left margin
#   ) +
#   # ggtitle(label = expression(paste(E*"["*n[PA]*"]", " = 800"))) +
#   ggtitle(label = expression(paste(E*"["*n[PO]*"]", " = 800 expected presence-only data"))) +
#   geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")
# 
# library(gridExtra)
# png(filename = paste0(getwd(), "/Figures/sim_results_spatially_constrained_lambda.png"), width = 4.2 * plot.res, height = 6.5 * plot.res, res = plot.res)
# grid.arrange(p1,p2,p3)
# dev.off()