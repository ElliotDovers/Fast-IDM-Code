## Analyse Scenario ##

home.wd <- getwd()

load("collated raw results.RDATA")

library(dplyr)
library(ggplot2)

all(paste0("PO", dat$scen_po, "_PA", dat$scen_pa) == dat$scen) # check I have assigned scenarios correctly
dat$scenario <- factor(paste0(dat$scen, "_", dat$env_range))

# we can look at the distribution of fitted likelihoods to understand which simulations result in poor model convergences
par(mfrow = c(2,2))
hist(dat$LL[dat$Model == "IDM" & !is.na(dat$LL)], main = "IDM", xlab = "logLik")
hist(dat$LL[dat$Model == "IDM" & dat$LL < 10000 & !is.na(dat$LL)], main = "IDM (logLik < 10,000)", xlab = "logLik")
hist(dat$LL[dat$Model == "PO Only" & !is.na(dat$LL)], main = "IDM", xlab = "logLik")
hist(dat$LL[dat$Model == "PA Only" & !is.na(dat$LL)], main = "IDM", xlab = "logLik")
par(mfrow = c(1,1))

# it appears that most of the positive log-Likelihoods are a poor convergence (occurring in the PO Only or IDM and none in the PA Only)
dat %>% filter(LL > 0 & !is.na(LL)) %>% group_by(scenario, Model) %>% summarise(KL_mu = sum(is.infinite(KLdiv)), KL_lambda = sum(is.infinite(KLdiv_lambda)))
# with just 3 instances of infinite KL Divergences otherwise
dat %>% filter(LL < 0 & !is.na(LL)) %>% group_by(scenario, Model) %>% summarise(KL_mu = sum(is.infinite(KLdiv)), KL_lambda = sum(is.infinite(KLdiv_lambda)))

# create full results table for Appendix
fit.issues <- dat %>% group_by(scen, env_range, latent_range, Model) %>%
  summarise(N_SIMS = length((sim)), INF_KL = sum(is.infinite(KLdiv)), INF_KL_lambda = sum(is.infinite(KLdiv_lambda)), POOR_CONV = sum(LL > 0 | is.nan(LL)))

# check that the infinite KLs are not concentrated so that there are not enough simulations
sim.check <- dat %>% filter(is.finite(KLdiv)) %>% group_by(Model, scen, env_range, bias_range) %>%
  summarise(unique_sims = length(unique(sim)), sims = length(sim), avg_k = mean(k), med_k = median(k), avg_k_bias = mean(k_bias), med_k_bias = median(k_bias))
min(sim.check$sims) # 84% is the lowest number returned in any one combination of env, lat and bias ranges
sim.check <- dat %>% filter(is.finite(KLdiv_lambda)) %>% group_by(Model, scen, env_range, bias_range) %>%
  summarise(unique_sims = length(unique(sim)), sims = length(sim), avg_k = mean(k), med_k = median(k), avg_k_bias = mean(k_bias), med_k_bias = median(k_bias))
min(sim.check$sims) # and 84% is also the case for KL divergence on the fitted lambda

# remove entire sims with infinite or NaN KLdiv (SO TO COMPARE MODELS FAIRLY)
dat$id <- paste(dat$scen, dat$env_range, dat$bias_range, dat$sim, sep = ":")
# first determine those sims for KL div on mu(s)
tmp <- dat %>% filter(Model != "PO Only") %>% group_by(scen, env_range, bias_range, sim) %>% summarise(check = sum(is.infinite(KLdiv) | is.na(KLdiv)))
rm.ids_mu <- paste(tmp$scen[tmp$check != 0], tmp$env_range[tmp$check != 0], tmp$bias_range[tmp$check != 0], tmp$sim[tmp$check != 0], sep = ":")
# then determine those sims for KL div on lambda(s)
tmp <- dat %>% filter(Model != "PA Only") %>% group_by(scen, env_range, bias_range, sim) %>% summarise(check = sum(is.infinite(KLdiv_lambda) | is.na(KLdiv_lambda)))
rm.ids_lambda <- paste(tmp$scen[tmp$check != 0], tmp$env_range[tmp$check != 0], tmp$bias_range[tmp$check != 0], tmp$sim[tmp$check != 0], sep = ":")
# these all occur in the same individual simulations
setequal(rm.ids_mu, rm.ids_lambda)
# remove the offending simulations for all models
dat <- dat[!dat$id %in% rm.ids_mu, ]

# there are also 4 KLdiv mu(s) results that are effectively Infinite that are throwing off results (3 for the PA Only model and 1 for IDM - which has an NaN logLik)
dat[dat$KLdiv > 1e7, ]
# which are similarly the only offending sims for KLdiv on lambda(s)
dat[dat$KLdiv > 1e7, ]
# going to remove them too
tmp <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(check = sum(KLdiv > 1e7))
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

# adjust our metric to be relative to the PA Only model for mu(s)
dat_pa <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(delta_kl = KLdiv[Model == "IDM"] - KLdiv[Model == "PA Only"])
dat_pa$scen_pa <- sapply(strsplit(dat_pa$scen, "_PA", fixed = T), function(x){x[2]})
dat_pa$Model <- paste0("IDM - N PO = ", sapply(strsplit(sapply(strsplit(dat_pa$scen, "PO", fixed = T), function(x){x[2]}), "_", fixed = T), function(x){x[1]}))
dat_pa$b_range <- factor(dat_pa$bias_range)
dat_pa$Model <- factor(dat_pa$Model, levels = c("IDM - N PO = 50", "IDM - N PO = 200", "IDM - N PO = 800"))
dat_pa$scen_pa <- factor(dat_pa$scen_pa, levels = c(50, 200, 800))

# adjust our metric to be relative to the PA Only model for lambda(s)
dat_po <- dat %>% group_by(scen, env_range, bias_range, sim) %>% summarise(delta_kl = KLdiv_lambda[Model == "IDM"] - KLdiv_lambda[Model == "PO Only"])
dat_po$scen_po <- sapply(strsplit(sapply(strsplit(dat_po$scen, "PO", fixed = T), function(x){x[2]}), "_", fixed = T), function(x){x[1]})
dat_po$Model <- paste0("IDM - N PA = ", sapply(strsplit(dat_po$scen, "_PA", fixed = T), function(x){x[2]}))
dat_po$b_range <- factor(dat_po$bias_range)
dat_po$Model <- factor(dat_po$Model, levels = c("IDM - N PA = 50", "IDM - N PA = 200", "IDM - N PA = 800"))
dat_po$scen_po <- factor(dat_po$scen_po, levels = c(50, 200, 800))

## Boxplots relative to the competing model ##

# set the ranges for each n[PA] and n[PO]
tmp <- boxplot(delta_kl ~ scen_pa, data = dat_pa)
pa_ranges <- data.frame(apply(tmp$stats, 2, function(x){c(floor(x[1]) + 0.05 * x[1], ceiling(x[5]) + 0.05 * x[5])}))
colnames(pa_ranges) <- c("pa50", "pa200", "pa800")
tmp <- boxplot(delta_kl ~ scen_po, data = dat_po)
po_ranges <- data.frame(apply(tmp$stats, 2, function(x){c(floor(x[1]) + 0.05 * x[1], ceiling(x[5]) + 0.05 * x[5])}))
colnames(po_ranges) <- c("po50", "po200", "po800")
# set some model colours
mod_cols_pa <- c("darkgoldenrod1", "darkgoldenrod3", "darkgoldenrod4")
mod_cols_po <- c("darkorchid2", "darkorchid3", "darkorchid4")

p1 <- ggplot(data = dat_pa %>% filter(scen_pa == 50 & b_range == "19"), aes(x = Model, y = delta_kl, color = Model, fill = Model)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected amount of PO integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")"))) +
  coord_cartesian(ylim = pa_ranges$pa50) +
  scale_color_manual(values = mod_cols_pa) +
  scale_fill_manual(values = alpha(mod_cols_pa, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", plot.title = element_text(size = 10),
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 12,  # Bottom margin
                             l = 0)  # Left margin
  ) +
  ggtitle(label = "Species expected at 5% of survey sites") +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p2 <- ggplot(data = dat_pa %>% filter(scen_pa == 200 & b_range == "19"), aes(x = Model, y = delta_kl, color = Model, fill = Model)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected amount of PO integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")"))) +
  coord_cartesian(ylim = pa_ranges$pa200) +
  scale_color_manual(values = mod_cols_pa) +
  scale_fill_manual(values = alpha(mod_cols_pa, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", plot.title = element_text(size = 10),
        plot.margin = margin(t = 12,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  ggtitle(label = "Species expected at 20% of survey sites") +
  # ggtitle(label = expression(paste("Species expected at 20%\nof survey sites (E["*n[PA]*"] = 200)"), sep = "")) +
  # ggtitle(label = expression(atop(Species~expected~at~'20%',of~survey~sites~(E*'['*n[PA]*']=200')))) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p3 <- ggplot(data = dat_pa %>% filter(scen_pa == 800 & b_range == "19"), aes(x = Model, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected amount of PO integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")), breaks = c(-20000, -10000, 0, 10000), labels = c("-20K", "-10K", "0", "10K")) +
  coord_cartesian(ylim = pa_ranges$pa800) +
  scale_color_manual(values = mod_cols_pa) +
  scale_fill_manual(values = alpha(mod_cols_pa, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=9), axis.title.x = element_text(margin = margin(t = 15)),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", plot.title = element_text(size = 10),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  ggtitle(label = "Species expected at 80% of survey sites") +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p4 <- ggplot(data = dat_po %>% filter(scen_po == 50 & b_range == "19"), aes(x = Model, y = delta_kl, color = Model, fill = Model)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected presences in PA integrated", labels = c(expression(paste(E*"["*n[PA]*"]", " = 50")), expression(paste(E*"["*n[PA]*"]", " = 200")), expression(paste(E*"["*n[PA]*"]", " = 800")))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")"))) +
  # coord_cartesian(ylim = po_ranges$po50) +
  coord_cartesian(ylim = c(-10,5)) +
  scale_color_manual(values = mod_cols_po) +
  scale_fill_manual(values = alpha(mod_cols_po, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 10),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = -1),
        plot.margin = margin(t = 20,  # Top margin
                             r = 0,  # Right margin
                             b = 12,  # Bottom margin
                             l = 0)  # Left margin
  ) +
  ggtitle(label = expression("50 expected presence-only records")) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p5 <- ggplot(data = dat_po %>% filter(scen_po == 200 & b_range == "19"), aes(x = Model, y = delta_kl, color = Model, fill = Model)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected presences in PA integrated", labels = c(expression(paste(E*"["*n[PO]*"]", " = 50")), expression(paste(E*"["*n[PO]*"]", " = 200")), expression(paste(E*"["*n[PO]*"]", " = 800")))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")"))) +
  coord_cartesian(ylim = po_ranges$po200) +
  scale_color_manual(values = mod_cols_po) +
  scale_fill_manual(values = alpha(mod_cols_po, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 10),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = -1),
        plot.margin = margin(t = 12,  # Top margin
                             r = 0,  # Right margin
                             b = 20,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  # ggtitle(label = expression(paste(E*"["*n[PO]*"]", " = 200 expected PO Data"))) +
  ggtitle(label = expression("200 expected presence-only records")) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

p6 <- ggplot(data = dat_po %>% filter(scen_po == 800 & b_range == "19"), aes(x = Model, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = "Expected presences in PA integrated", labels = c(expression(paste(E*"["*n[PA]*"]", " = 50")), expression(paste(E*"["*n[PA]*"]", " = 200")), expression(paste(E*"["*n[PA]*"]", " = 800")))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")"))) +
  coord_cartesian(ylim = po_ranges$po800) +
  scale_color_manual(values = mod_cols_po) +
  scale_fill_manual(values = alpha(mod_cols_po, alpha = 0.5)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 10),
        panel.background = element_rect(fill = NA, color = "black"), axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=10), axis.text.x = element_text(size=9), axis.title.x = element_text(margin = margin(t = 15)),
        legend.text=element_text(size=10), legend.text.align = 0, legend.position = "none", axis.title.y = element_text(vjust = -1),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)  # Left margin
  ) + 
  ggtitle(label = expression("800 expected presence-only records")) +
  geom_abline(intercept = 0, slope = 0, lty = "dashed", color = "darkblue")

library(gridExtra)
plot.res <- 500
png(filename = paste0(home.wd, "/Figures/our_sim_results_boxplots.png"), width = 6.2 * plot.res, height = 5.5 * plot.res, res = plot.res)
grid.arrange(p1,p4,p2,p5,p3,p6, nrow = 3, ncol = 2)
dev.off()

################################################################################
### Appendix plot - boxplots across bias ranges ################################

# set a scale function to round to KL div. axis values to zero decimal places
scaleFUN <- function(x) sprintf("%.0f", x)

Ap1 <- ggplot(data = dat_pa %>% filter(scen_pa == 50 & env_range == 15), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s)))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = c(-100, 50)) +
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
  coord_cartesian(ylim = c(-120, 70)) +
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
  coord_cartesian(ylim = c(-20000, 11000)) +
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
  coord_cartesian(ylim = c(-10, 5)) +
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
  coord_cartesian(ylim = c(-10, 5)) +
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
  coord_cartesian(ylim = c(-12, 5)) +
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
png(filename = paste0(home.wd, "/Figures/our_sim_results_boxplots_append_short_env.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
grid.arrange(Ap1,Ap4,Ap2,Ap5,Ap3,Ap6, nrow = 3, ncol = 2)
dev.off()

Ap1 <- ggplot(data = dat_pa %>% filter(scen_pa == 50 & env_range == 25), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s)))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = c(-110, 50)) +
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

Ap2 <- ggplot(data = dat_pa %>% filter(scen_pa == 200 & env_range == 25), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s)))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = c(-220, 90)) +
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

Ap3 <- ggplot(data = dat_pa %>% filter(scen_pa == 800 & env_range == 25), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s))) , breaks = seq(11, 29, by = 4)) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", mu, " || ", hat(mu), ")")), breaks = c(-30000, -20000, -10000, 0, 10000), labels = c("-30K", "-20K", "-10K", "0", "10K")) +
  coord_cartesian(ylim = c(-38000, 16000)) +
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

Ap4 <- ggplot(data = dat_po %>% filter(scen_po == 50 & env_range == 25), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s)))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = c(-12, 5)) +
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

Ap5 <- ggplot(data = dat_po %>% filter(scen_po == 200 & env_range == 25), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s)))) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = c(-11, 5)) +
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

Ap6 <- ggplot(data = dat_po %>% filter(scen_po == 800 & env_range == 25), aes(x = b_range, y = delta_kl, col = Model, fill = Model)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_boxplot(outlier.shape = 1, outlier.size = 1) +
  scale_x_discrete(name = expression(paste(rho[B], ": Spatial range parameter of ", xi[B](s))) , breaks = seq(11, 29, by = 4)) +
  scale_y_continuous(name = expression(paste(Delta~D[KL],"(", lambda, " || ", hat(lambda), ")")), labels=scaleFUN) +
  coord_cartesian(ylim = c(-13, 5)) +
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
png(filename = paste0(home.wd, "/Figures/our_sim_results_boxplots_append_long_env.png"), width = 6.2 * plot.res, height = 6.2 * plot.res, res = plot.res)
grid.arrange(Ap1,Ap4,Ap2,Ap5,Ap3,Ap6, nrow = 3, ncol = 2)
dev.off()

################################################################################
# Append Full Results Table ####################################################
################################################################################

tab <- dat %>% group_by(scen, env_range, latent_range, Model) %>%
  summarise(KL = median(KLdiv), KL_lambda = median(KLdiv_lambda), Comp_Time = median(timing))
# set up a unique id to link to fit issues table
id1 <- paste(tab$scen, tab$env_range, tab$latent_range, tab$Model, sep = ":")
id2 <- paste(fit.issues$scen, fit.issues$env_range, fit.issues$latent_range, fit.issues$Model, sep = ":")
# add in the fit issues
tab <- cbind(tab, fit.issues[match(id1, id2), c("N_SIMS", "INF_KL", "INF_KL_lambda", "POOR_CONV")])

# separate the expected sample size scenarios
tab$n_pa <- factor(sapply(strsplit(tab$scen, "PA", fixed = T), function(x){x[2]}), levels = c("50", "200", "800"))
tab$n_po <- factor(sapply(strsplit(sapply(strsplit(tab$scen, "_", fixed = T), function(x){x[1]}), "PO", fixed = T), function(x){x[2]}), levels = c("50", "200", "800"))

# re-order and align the columns for presenting
tab <- tab[order(tab$n_pa, tab$n_po), c("n_pa", "n_po", "env_range", "latent_range", "Model", "KL", "KL_lambda", "Comp_Time", "INF_KL", "INF_KL_lambda", "POOR_CONV")]

library(xtable)
print(xtable(tab, label = "append:tab:sim", digits = 2,
             caption = c("Full results from our simulations")),
      include.rownames=F, caption.placement = "top")


################################################################################