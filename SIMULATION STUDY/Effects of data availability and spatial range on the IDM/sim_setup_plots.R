home.wd <- getwd()
# Get the job array
tab <- read.csv("job_array.csv")

# determine job number from pbs script
job = 701

################################################################################
# Set parameters that define the scenarios:

# random seed / simulation number
seed = tab$sim_number[tab$job == job]
# environment range of effect
env.range = tab$env_corrleation_range[tab$job == job]
# latent range of effect
lat.range = tab$lat_corrleation_range[tab$job == job]
# bias field range of effect
bias.range = tab$bias_corrleation_range[tab$job == job]

################################################################################

source("sim_PO_PA_data.R")
dat_pa <- sim_occurrence_data(Intercept_po = -4,
                                       Intercept_pa = -1.25,
                                       rseed = seed,
                                       env.covariate.type = "random_field",
                                       presence.only.observer.bias.covariate.type = "random_field",
                                       presence.only.observer.bias.covariate.range = bias.range,
                                       env.covariate.range = env.range,
                                       latent.range = lat.range,
                                       latent.field = T,
                                       plotting = T
                                       
)
dat_po <- attr(dat_pa, "presence-only")
quad <- attr(dat_pa, "truth.grid")

dat_pa2 <- sim_occurrence_data(Intercept_po = -4,
                              Intercept_pa = -1.25,
                              rseed = 10,
                              env.covariate.type = "random_field",
                              presence.only.observer.bias.covariate.type = "random_field",
                              presence.only.observer.bias.covariate.range = 15,
                              env.covariate.range = env.range,
                              latent.range = lat.range,
                              latent.field = T,
                              plotting = T
                              
)
dat_po2 <- attr(dat_pa2, "presence-only")
quad2 <- attr(dat_pa2, "truth.grid")

plot.res <- 500

png(filename = paste0(getwd(), "/Results/Figures/our_sim_setup.png"), res = plot.res, width = 5 * plot.res, height = 5 * plot.res)
layout(mat = matrix(c(1,2,3,4,5,6,7,8,9), nrow = 3, ncol = 3, byrow = F), heights = rep(1/3,3), widths = rep(1/3,3))

par(mar = rep(3, 4))
# Env. Covariate
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(quad$env, quad$x, quad$y),
      col = terrain.colors(256), xlab = "", ylab = "", main = expression(paste("Measured Effects: ", bold(X))), bty = 'n',
      axes = F, asp = 1
)
lines(c(1, 1, 100, 100, 1), c(1, 100, 100, 1, 1))

# Latent Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(quad$xi[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$x[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$y[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]),
      col = plasma(256), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
title(main = expression(paste("Unmeasured Effects: ", xi)), line = 1)
symbols(x = 50.5, y = 50.5, circles = 50.5, inches = F, add = T)

# Dummy
plot(1, type = "n", axes = F)

# Abundance Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(log(quad$mu[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]), quad$x[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$y[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]),
      col = rainbow(256), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
title(main = expression(paste("Abundance Rate: ", mu)), line = 1)
symbols(x = 50.5, y = 50.5, circles = 50.5, inches = F, add = T)


# Bias Covariate
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(quad$bias[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$x[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$y[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]),
      col = cividis(256), xlab = "", ylab = "", main = expression(paste("PO Biasing: ", B)), bty = 'n',
      axes = F, asp = 1
)
symbols(x = 50.5, y = 50.5, circles = 50.5, inches = F, add = T)

# Intensity Field
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(log(quad$lambda[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]), quad$x[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$y[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]),
      col = rainbow(256), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
title(main = expression(paste("Pres. Record Rate: ", lambda)), line = 1)
symbols(x = 50.5, y = 50.5, circles = 50.5, inches = F, add = T)


# PA data
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(log(quad$mu), quad$x, quad$y),
      col = NA, xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
title(main = expression(paste(PA, " Data")), line = 1)
pa.col <- dat_pa$pa
pa.col[dat_pa$present == 0] <- "lightblue"
pa.col[dat_pa$present == 1] <- "darkblue"
pa.pch <- dat_pa$present
pa.pch[dat_pa$present == 0] <- 4
pa.pch[dat_pa$present == 1] <- 1
points(dat_pa[ , c("x", "y")], pch = pa.pch, col = pa.col)
lines(c(1, 1, 100, 100, 1), c(1, 100, 100, 1, 1))

# Dummy
plot(1, type = "n", axes = F)
legend("center", legend = c("PA: present", "PA: absence", "PO: presence"),
       pch = c(1, 4, 1), col = c("darkblue", "lightblue", "darkred"), bty = "n"
)

# PO data
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(log(quad$lambda), quad$x, quad$y),
      col = NA, xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
title(main = expression(paste(PO, " Data")), line = 1)
points(dat_po[ , c("x", "y")], col = "darkred")
lines(c(1, 1, 100, 100, 1), c(1, 100, 100, 1, 1))

dev.off()

### Individual plots

# Env. Covariate
png(filename = paste0(getwd(), "/Results/Figures/sim_setup_X.png"), res = plot.res, width = 5 * plot.res, height = 5 * plot.res)
par(mar = rep(0,4), oma = rep(0,4))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(quad$env, quad$x, quad$y),
      col = terrain.colors(256, alpha = 0.6), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
lines(c(1, 1, 100, 100, 1), c(1, 100, 100, 1, 1))
text(50, 50, labels = "X(s)", cex = 10, adj = 0.5)
dev.off()

# Latent Field
png(filename = paste0(getwd(), "/Results/Figures/sim_setup_xi.png"), res = plot.res, width = 5 * plot.res, height = 5 * plot.res)
par(mar = rep(0,4), oma = rep(0,4))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(quad$xi[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$x[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$y[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]),
      col = plasma(256, alpha = 0.6), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
symbols(x = 50.5, y = 50.5, circles = 50.5, inches = F, add = T)
text(50, 50, labels = expression(bold(xi)(s)), cex = 10, adj = 0.5)
dev.off()


# Abundance Field
png(filename = paste0(getwd(), "/Results/Figures/sim_setup_abund.png"), res = plot.res, width = 5 * plot.res, height = 5 * plot.res)
par(mar = rep(0,4), oma = rep(0,4))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(log(quad$mu[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]), quad$x[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$y[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]),
      col = rainbow(256, alpha = 0.6), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
symbols(x = 50.5, y = 50.5, circles = 50.5, inches = F, add = T)
text(50, 50, labels = expression(bold(mu)(s)), cex = 10, adj = 0.5)
dev.off()


# Bias Covariate
png(filename = paste0(getwd(), "/Results/Figures/sim_setup_bias.png"), res = plot.res, width = 5 * plot.res, height = 5 * plot.res)
par(mar = rep(0,4), oma = rep(0,4))
image(x = sort(unique(quad2$x)),
      y = sort(unique(quad2$y)),
      z = scampr:::vec2mat(quad2$bias, quad2$x, quad2$y),
      col = cividis(256, alpha = 0.6), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
lines(c(1, 1, 100, 100, 1), c(1, 100, 100, 1, 1))
text(50, 50, labels = expression(bold(X)[B](s)), cex = 10, adj = 0.5)
dev.off()

# Bias Latent Field
png(filename = paste0(getwd(), "/Results/Figures/sim_setup_bias_xi.png"), res = plot.res, width = 5 * plot.res, height = 5 * plot.res)
par(mar = rep(0,4), oma = rep(0,4))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(quad$bias[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$x[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$y[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]),
      col = cividis(256, alpha = 0.6), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
symbols(x = 50.5, y = 50.5, circles = 50.5, inches = F, add = T)
text(50, 50, labels = expression(bold(xi)[B](s)), cex = 10, adj = 0.5)
dev.off()

# Intensity Field
png(filename = paste0(getwd(), "/Results/Figures/sim_setup_inten.png"), res = plot.res, width = 5 * plot.res, height = 5 * plot.res)
par(mar = rep(0,4), oma = rep(0,4))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(log(quad$lambda[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]), quad$x[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2], quad$y[(quad$x-50.5)^2 + (quad$y-50.5)^2 <= 50.5^2]),
      col = rainbow(256, alpha = 0.6), xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
symbols(x = 50.5, y = 50.5, circles = 50.5, inches = F, add = T)
text(50, 50, labels = expression(bold(lambda)(s)), cex = 10, adj = 0.5)
dev.off()


# PA data
png(filename = paste0(getwd(), "/Results/Figures/sim_setup_pa.png"), res = plot.res, width = 5 * plot.res, height = 5 * plot.res)
par(mar = rep(0,4), oma = rep(0,4))
# par(mar = c(0,0,2.1,0), oma = rep(0,4))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(log(quad$mu), quad$x, quad$y),
      col = NA, xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
pa.col <- dat_pa$pa
pa.col[dat_pa$present == 0] <- "darkblue"
pa.col[dat_pa$present == 1] <- "darkred"
pa.pch <- dat_pa$present
pa.pch[dat_pa$present == 0] <- 4
pa.pch[dat_pa$present == 1] <- 16
points(dat_pa[ , c("x", "y")], pch = pa.pch, col = pa.col, cex = 1.5)
lines(c(1, 1, 100, 100, 1), c(1, 100, 100, 1, 1))
text(50, 50, labels = "Presence\nAbsence\nData", cex = 6, adj = 0.5)
dev.off()

# PO data
png(filename = paste0(getwd(), "/Results/Figures/sim_setup_po.png"), res = plot.res, width = 5 * plot.res, height = 5 * plot.res)
par(mar = rep(0,4), oma = rep(0,4))
image(x = sort(unique(quad$x)),
      y = sort(unique(quad$y)),
      z = scampr:::vec2mat(log(quad$lambda), quad$x, quad$y),
      col = NA, xlab = "", ylab = "", main = "", bty = 'n',
      axes = F, asp = 1
)
points(dat_po[ , c("x", "y")], col = "darkred", pch = 1, cex = 2)
lines(c(1, 1, 100, 100, 1), c(1, 100, 100, 1, 1))
text(50, 50, labels = "Presence\nOnly\nData", cex = 6, adj = 0.5)
dev.off()